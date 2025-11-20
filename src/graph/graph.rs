use hashbrown::HashMap;
use petgraph::{
    algo::maximal_cliques::maximal_cliques,
    graph::{Graph, NodeIndex},
    Undirected,
};

use duckdb::{Connection, Result, params};

/// Returns `Some((nodes, edges))` if there is a clique ≥ `click_size`, or `None` if not.
/// * `signatures` – list of signatures (u64) to consider.
/// * `dist_max` – maximum threshold for the `dist` column.
/// * `click_size` – minimum size of the clique to be found.
pub fn exists_clique(
    conn: &Connection,
    signatures: &[u64],
    dist_max: f64,
    click_size: usize,
) -> Result<Option<(Vec<u64>, Vec<u64>)>> {
    if signatures.is_empty() || click_size == 0 {
        return Ok(None);
    }

    /* ---------- 1. Filtered SQL query ---------- */
    let placeholders = signatures
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

    /* parameters: signatures for source, signatures for target, dist_max */
    let mut params: Vec<&dyn duckdb::ToSql> = Vec::new();
    
    // Add signatures for source IN clause
    for sig in signatures {
        params.push(sig);
    }
    // Add signatures for target IN clause
    for sig in signatures {
        params.push(sig);
    }
    params.push(&dist_max);

    let mut stmt = conn.prepare(&sql)?;
    let rows = stmt.query_map(params.as_slice(), |row| {
        Ok((
            row.get::<_, u64>(0)?,     // row id
            row.get::<_, u64>(1)?,     // source (BIGINT)
            row.get::<_, u64>(2)?,     // target (BIGINT)
        ))
    })?;

    /* ---------- 2. Graph construction ---------- */
    let mut g: Graph<u64, u64, Undirected> = Graph::new_undirected();
    let mut sig2idx: HashMap<u64, NodeIndex> = HashMap::new();

    for &sig in signatures {
        let n = g.add_node(sig);
        sig2idx.insert(sig, n);
    }

    let mut edge_lookup: HashMap<(u64, u64), u64> = HashMap::new();

    for row in rows.flatten() {
        let (edge_id, s, t) = row;
        let s_u64 = s as u64;
        let t_u64 = t as u64;
        
        if let (Some(&a), Some(&b)) = (sig2idx.get(&s_u64), sig2idx.get(&t_u64)) {
            g.update_edge(a, b, edge_id);
            let key = if s_u64 <= t_u64 { 
                (s_u64, t_u64) 
            } else { 
                (t_u64, s_u64) 
            };
            edge_lookup.insert(key, edge_id);
        }
    }

    /* ---------- 3. Clique search ---------- */
    for clique in maximal_cliques(&g) {
        if clique.len() >= click_size {
            /* Retrieve node signatures */
            let mut nodes: Vec<u64> = clique.iter().map(|&ix| g[ix]).collect();
            nodes.sort();

            /* Retrieve edge IDs from the clique */
            let mut edges = Vec::<u64>::new();
            for i in 0..nodes.len() {
                for j in (i + 1)..nodes.len() {
                    let key = if nodes[i] <= nodes[j] {
                        (nodes[i], nodes[j])
                    } else {
                        (nodes[j], nodes[i])
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
/// - `edge_ids`: Slice of `u64` IDs (primary key) to delete.
/// - `dist_min`: Minimum distance; only edges with `dist < dist_min` will be deleted.
///
/// # Returns:
/// - `Ok(())` on success, or `Err` on database error.
pub fn delete_edges_by_ids(
    conn: &mut Connection,
    edge_ids: &[u64],
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

/// Gets edges from `edges` where source or target is in the signatures slice (both BIGINT)
/// and dist < dist_max. Returns a vector of (source, target, dist).
pub fn get_edges_by_signatures_and_dist(
    conn: &Connection,
    signatures: &[u64],
    dist_max: f64,
) -> Result<Vec<(u64, u64, f64)>> {
    if signatures.is_empty() {
        return Ok(Vec::new());
    }
    
    // placeholders ?,?,?,...
    let placeholders = signatures.iter().map(|_| "?").collect::<Vec<_>>().join(", ");
    let sql = format!(
        "SELECT source, target, dist FROM edges \
        WHERE (source IN ({}) OR target IN ({})) AND dist < ?",
        placeholders, placeholders
    );

    // Parameters: signatures for source, signatures for target, dist_max
    let mut params: Vec<&dyn duckdb::ToSql> = Vec::new();
    
    // Add signatures for source IN clause
    for sig in signatures {
        params.push(sig);
    }
    // Add signatures for target IN clause
    for sig in signatures {
        params.push(sig);
    }
    params.push(&dist_max);

    let mut stmt = conn.prepare(&sql)?;
    let mut result = Vec::new();

    let rows = stmt.query_map(params.as_slice(), |row| {
        Ok((
            row.get::<_, u64>(0)? as u64, // source (BIGINT -> u64)
            row.get::<_, u64>(1)? as u64, // target (BIGINT -> u64)
            row.get::<_, f64>(2)?          // dist
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
/// - edges: slice of tuples (source signature, target signature, dist)
/// 
/// # Note
/// - The table must have columns source (BIGINT), target (BIGINT), dist (DOUBLE)
pub fn insert_edges(
    conn: &Connection,
    edges: &[(u64, u64, f64)]
) -> Result<()> {
    let sql = "INSERT INTO edges (source, target, dist) VALUES (?, ?, ?)";
    let mut stmt = conn.prepare(sql)?;

    for (source, target, dist) in edges {
        stmt.execute(params![*source as u64, *target as u64, dist])?;
    }
    Ok(())
}
