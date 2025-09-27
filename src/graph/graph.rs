use hashbrown::HashMap;
use petgraph::{
    algo::maximal_cliques::maximal_cliques,
    graph::{Graph, NodeIndex},
    Undirected,
};

use duckdb::{Connection, Result, params};
/// Returns `Some((nodes, edges))` if there is a clique ≥ `click_size`, or `None` if not.
/// * `ids`      – list of IDs (VARCHAR) to consider.
/// * `dist_max` – maximum threshold for the `dist` column.
/// * `click_size` – minimum size of the clique to be found.
pub fn exists_clique(
    conn: &Connection,
    ids: &[String],
    dist_max: f64,
    click_size: usize,
) -> Result<Option<(Vec<String>, Vec<i64>)>> {
    if ids.is_empty() || click_size == 0 {
        return Ok(None);
    }

    /* ---------- 1. Filtered SQL query ---------- */
    let placeholders = ids
        .iter()
        .map(|_| "?")
        .collect::<Vec<_>>()
        .join(", ");
    let sql = format!(
        "SELECT id, source, target       \
         FROM edges                     \
         WHERE source IN ({})           \
           AND target IN ({})           \
           AND dist > ?",
        placeholders, placeholders
    );

    /* parameters: ids for source, ids for target, dist_max */
    let mut params: Vec<&dyn duckdb::ToSql> =
        ids.iter().map(|s| s as &dyn duckdb::ToSql).collect();
    params.extend(ids.iter().map(|s| s as &dyn duckdb::ToSql));
    params.push(&dist_max);

    let mut stmt = conn.prepare(&sql)?;
    let rows = stmt.query_map(params.as_slice(), |row| {
        Ok((
            row.get::<_, i64>(0)?,     // row id
            row.get::<_, String>(1)?,  // source
            row.get::<_, String>(2)?,  // target
        ))
    })?;

    /* ---------- 2. Graph construction ---------- */
    let mut g: Graph<String, i64, Undirected> = Graph::new_undirected();
    let mut id2idx: HashMap<String, NodeIndex> = HashMap::new();

    for id in ids {
        let n = g.add_node(id.clone());
        id2idx.insert(id.clone(), n);
    }

    let mut edge_lookup: HashMap<(String, String), i64> = HashMap::new();

    for row in rows.flatten() {
        let (edge_id, s, t) = row;
        if let (Some(&a), Some(&b)) = (id2idx.get(&s), id2idx.get(&t)) {
            g.update_edge(a, b, edge_id);
            let key = if s <= t { (s.clone(), t) } else { (t, s) };
            edge_lookup.insert(key, edge_id);
        }
    }

    /* ---------- 3. Clique search ---------- */
    for clique in maximal_cliques(&g) {
        if clique.len() >= click_size {
            /* Retrieve node IDs */
            let mut nodes: Vec<String> = clique.iter().map(|&ix| g[ix].clone()).collect();
            nodes.sort();

            /* Retrieve edge IDs from the clique */
            let mut edges = Vec::<i64>::new();
            for i in 0..nodes.len() {
                for j in (i + 1)..nodes.len() {
                    let key = if nodes[i] <= nodes[j] {
                        (nodes[i].clone(), nodes[j].clone())
                    } else {
                        (nodes[j].clone(), nodes[i].clone())
                    };
                    if let Some(&eid) = edge_lookup.get(&key) {
                        edges.push(eid);
                    }
                }
            }
            return Ok(Some((nodes, edges)));
        }
    }

    Ok(None)
}

/// Deletes from the 'edges' table all rows whose ID is in `edge_ids`
/// and whose distance is less than `dist_min`.
///
/// # Parameters:
/// - `conn`: mutable DuckDB connection.
/// - `edge_ids`: Slice of `i64` IDs (primary key) to delete.
/// - `dist_min`: Minimum distance; only edges with `dist < dist_min` will be deleted.
///
/// # Returns:
/// - `Ok(())` on success, or `Err` on database error.
pub fn delete_edges_by_ids(
    conn: &mut Connection,
    edge_ids: &[i64],
    dist_min: f64,
) -> Result<(), duckdb::Error> {
    if edge_ids.is_empty() {
        return Ok(());
    }
    
    let tx = conn.transaction()?;
    
    // Create a temporary table with a unique name
    let temp_table = format!("temp_ids_{}", std::process::id());
    tx.execute_batch(&format!("CREATE TEMP TABLE {} (id INTEGER)", temp_table))?;
    
    // Insert IDs into the temporary table
    let mut insert_stmt = tx.prepare(&format!("INSERT INTO {} (id) VALUES (?)", temp_table))?;
    for &id in edge_ids {
        insert_stmt.execute(params![id])?;
    }
    
    // Delete only those edges whose id is in temp_table and dist < dist_min
    let deleted = tx.execute(
        &format!(
            "DELETE FROM edges \
             WHERE id IN (SELECT id FROM {}) \
               AND dist < ?",
            temp_table
        ),
        params![dist_min],
    )?;
    
    // Clean up the temporary table
    tx.execute(&format!("DROP TABLE {}", temp_table), [])?;
    
    tx.commit()?;
    
    println!("✓ Deleted {} rows from edges table", deleted);
    Ok(())
}

/// Gets edges from `edges` where source or target is in the node_ids slice (both VARCHAR)
/// and dist < dist_max. Returns a vector of (source, target, dist).
pub fn get_edges_by_node_ids_and_dist(
    conn: &Connection,
    node_ids: &[&str],
    dist_max: f64,
) -> Result<Vec<(String, String, f64)>> {
    if node_ids.is_empty() {
        return Ok(Vec::new());
    }
    // placeholders ?,?,?,...
    let placeholders = node_ids.iter().map(|_| "?").collect::<Vec<_>>().join(", ");
    let sql = format!(
        "SELECT source, target, dist FROM edges \
        WHERE (source IN ({}) OR target IN ({})) AND dist < ?",
        placeholders, placeholders
    );

    // Parameters: node_ids for source, node_ids for target, dist_max
    let mut params: Vec<&dyn duckdb::ToSql> = node_ids.iter().map(|id| id as &dyn duckdb::ToSql).collect();
    params.extend(node_ids.iter().map(|id| id as &dyn duckdb::ToSql));
    params.push(&dist_max);

    let mut stmt = conn.prepare(&sql)?;
    let mut result = Vec::new();

    let rows = stmt.query_map(params.as_slice(), |row| {
        Ok((
            row.get::<_, String>(0)?, // source
            row.get::<_, String>(1)?, // target
            row.get::<_, f64>(2)?     // dist
        ))
    })?;

    for row in rows.flatten() {
        result.push(row);
    }
    Ok(result)
}

/// Inserts a vector of edges (source, target, dist) into the DuckDB edges table.
/// 
/// # Parameters
/// - conn: active DuckDB connection
/// - edges: slice of tuples (source, target, dist)
/// 
/// # Note
/// - The table must have columns source (VARCHAR), target (VARCHAR), dist (DOUBLE)
pub fn insert_edges(
    conn: &Connection,
    edges: &[(String, String, f64)]
    ) -> Result<()> {
    let sql = "INSERT INTO edges (source, target, dist) VALUES (?, ?, ?)";
    let mut stmt = conn.prepare(sql)?;

    for (source, target, dist) in edges {
        stmt.execute(params![source, target, dist])?;
    }
    Ok(())
}
