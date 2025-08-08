use hashbrown::HashMap;
use petgraph::{
    algo::maximal_cliques::maximal_cliques,
    graph::{Graph, NodeIndex},
    Undirected,
};
use duckdb::{Connection, Result, params};

/// Devuelve `Some((nodos, edges))` si hay un clique ≥ `click_size`, o `None` si no lo hay.
/// * `ids`      – lista de IDs (VARCHAR) a considerar.
/// * `dist_max` – umbral máximo para la columna `dist`.
/// * `click_size` – tamaño mínimo del clique buscado.
pub fn exists_clique(
    conn: &Connection,
    ids: &[String],
    dist_max: f64,
    click_size: usize,
) -> Result<Option<(Vec<String>, Vec<i64>)>> {
    if ids.is_empty() || click_size == 0 {
        return Ok(None);
    }

    /* ---------- 1. Consulta SQL filtrada ---------- */
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

    /* parámetros: ids para source, ids para target, dist_max */
    let mut params: Vec<&dyn duckdb::ToSql> =
        ids.iter().map(|s| s as &dyn duckdb::ToSql).collect();
    params.extend(ids.iter().map(|s| s as &dyn duckdb::ToSql));
    params.push(&dist_max);

    let mut stmt = conn.prepare(&sql)?;
    let rows = stmt.query_map(params.as_slice(), |row| {
        Ok((
            row.get::<_, i64>(0)?,     // id fila
            row.get::<_, String>(1)?,  // source
            row.get::<_, String>(2)?,  // target
        ))
    })?;

    /* ---------- 2. Construcción del grafo ---------- */
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

    /* ---------- 3. Búsqueda de cliques ---------- */
    for clique in maximal_cliques(&g) {
        if clique.len() >= click_size {
            /* Recuperar IDs de nodos */
            let mut nodes: Vec<String> = clique.iter().map(|&ix| g[ix].clone()).collect();
            nodes.sort();

            /* Recuperar IDs de edges del clique */
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


/// Elimina de la tabla 'edges' todas las filas cuyo ID esté en `edge_ids`
/// y cuya distancia sea menor que `dist_max`.
///
/// # Parámetros:
/// - `conn`      : Conexión mutable a DuckDB.
/// - `edge_ids`  : Slice de `i64` con los IDs (primary key) a eliminar.
/// - `dist_max`  : Distancia minima; solo se eliminan edges con `dist < dist_min`.
///
/// # Retorna:
/// - `Ok(())` en caso de éxito, o `Err` en caso de error de base de datos.
pub fn delete_edges_by_ids(
    conn: &mut Connection,
    edge_ids: &[i64],
    dist_min: f64,
) -> Result<(), duckdb::Error> {
    if edge_ids.is_empty() {
        return Ok(());
    }
    
    let tx = conn.transaction()?;
    
    // Crear tabla temporal con nombre único
    let temp_table = format!("temp_ids_{}", std::process::id());
    tx.execute_batch(&format!("CREATE TEMP TABLE {} (id INTEGER)", temp_table))?;
    
    // Insertar IDs en la tabla temporal
    let mut insert_stmt = tx.prepare(&format!("INSERT INTO {} (id) VALUES (?)", temp_table))?;
    for &id in edge_ids {
        insert_stmt.execute(params![id])?;
    }
    
    // Eliminar solo aquellos edges con id en temp_table y dist < dist_min
    let deleted = tx.execute(
        &format!(
            "DELETE FROM edges \
             WHERE id IN (SELECT id FROM {}) \
               AND dist < ?",
            temp_table
        ),
        params![dist_min],
    )?;
    
    // Limpiar tabla temporal
    tx.execute(&format!("DROP TABLE {}", temp_table), [])?;
    
    tx.commit()?;
    
    println!("✓ Eliminadas {} filas de la tabla edges", deleted);
    Ok(())
}


/// Obtiene los edges de `edges` donde source o target está en el slice de node_ids (ambos VARCHAR)
/// y dist < dist_max. Retorna un vector de (source, target, dist).
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

    // Parametros: node_ids para source, node_ids para target, dist_max
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


/// Inserta un vector de edges (source, target, dist) en la tabla edges de DuckDB.
/// 
/// # Parámetros
/// - conn: conexión activa a DuckDB
/// - edges: slice de tuplas (source, target, dist)
/// 
/// # Nota
/// - La tabla debe tener columnas source (VARCHAR), target (VARCHAR), dist (DOUBLE)
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

