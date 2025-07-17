use rayon::prelude::*;
use hashbrown::HashMap;
use petgraph::graph::{Graph, NodeIndex};
use petgraph::Undirected;
use petgraph::algo::clique::BronKerboschAllCliques;
use std::collections::HashMap;
use duckdb::{Connection, Result, params};

/// Busca un clique y devuelve los nodos y los IDs únicos de los edges usados
///
/// # Parámetros:
/// - conn: &Connection (DuckDB)
/// - ids: &[i64]        -- IDs a considerar como nodos
/// - dist_max: f64      -- valor máximo de Dist
/// - click_size: usize  -- tamaño mínimo del clique buscado
///
/// # Devuelve:
/// - Option<(Vec<i64>, Vec<i64>)> : (IDs de nodos, IDs de edges usados), o None si no hay clique
pub fn exists_clique(
    conn: &Connection,
    ids: &[i64],
    dist_max: f64,
    click_size: usize
) -> Result<Option<(Vec<i64>, Vec<i64>)>> {
    if ids.is_empty() || click_size == 0 {
        return Ok(None);
    }

    // 1. Prepara IDs para consulta SQL
    let id_list = ids.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(",");
    let query = format!(
        "SELECT id, source, target FROM edges \
         WHERE source IN ({0}) AND target IN ({0}) AND dist < ?",
         id_list
    );
    let mut stmt = conn.prepare(&query)?;
    let rows = stmt.query_map([dist_max], |row| {
        let row_id: i64 = row.get(0)?;
        let src: i64 = row.get(1)?;
        let tgt: i64 = row.get(2)?;
        Ok((row_id, src, tgt))
    })?;

    // 2. Arma el grafo y un lookup de aristas -> IDs originales
    let mut id_to_index: HashMap<i64, NodeIndex> = HashMap::new();
    let mut index_to_id: Vec<i64> = Vec::new();
    let mut graph = Graph::<i64, i64, Undirected>::new_undirected();
    // Nodos
    for &id in ids {
        let idx = graph.add_node(id);
        id_to_index.insert(id, idx);
        index_to_id.push(id);
    }
    // Mapea (src, tgt) -> edge_id
    let mut edge_id_lookup: HashMap<(i64, i64), i64> = HashMap::new();
    // Aristas
    for row in rows.flatten() {
        let (row_id, src, tgt) = row;
        if let (Some(&a), Some(&b)) = (id_to_index.get(&src), id_to_index.get(&tgt)) {
            graph.update_edge(a, b, row_id);
            // Para consulta rápida posterior
            let mut key = (src, tgt);
            // Siempre ordena el par (por si el grafo es no dirigido)
            if key.0 > key.1 { std::mem::swap(&mut key.0, &mut key.1); }
            edge_id_lookup.insert(key, row_id);
        }
    }
    // 3. Busca los cliques
    for clique_indices in BronKerboschAllCliques::new(&graph).filter(|c| c.len() >= click_size) {
        // IDs de los nodos del clique
        let clique_node_ids: Vec<i64> = clique_indices.iter().map(|ix| graph[*ix]).collect();
        // IDs de los edges usados en ese clique (todas las conexiones pares dentro del clique)
        let mut clique_edge_ids = Vec::new();
        for i in 0..clique_node_ids.len() {
            for j in (i+1)..clique_node_ids.len() {
                let mut key = (clique_node_ids[i], clique_node_ids[j]);
                if key.0 > key.1 { std::mem::swap(&mut key.0, &mut key.1); }
                if let Some(&eid) = edge_id_lookup.get(&key) {
                    clique_edge_ids.push(eid);
                }
            }
        }
        return Ok(Some((clique_node_ids, clique_edge_ids)));
    }
    Ok(None)
}



/// Elimina de la tabla 'edges' todas las filas cuyo ID esté en `edge_ids`.
///
/// # Parámetros:
/// - conn: Conexión DuckDB.
/// - edge_ids: Slice de i64 con los IDs (primary key) a eliminar.
///
/// # Retorna:
/// - Ok(()) en caso de éxito, o Err en caso de error de base de datos.
///
/// # Nota:
/// Si el vector está vacío, la función no hace nada.
pub fn delete_edges_by_ids(conn: &Connection, edge_ids: &[i64]) -> Result<()> {
    if edge_ids.is_empty() {
        return Ok(());
    }
    // Construye la lista de "?, ?, ?, ..." según el número de elementos
    let placeholders = edge_ids.iter().map(|_| "?").collect::<Vec<_>>().join(", ");
    let sql = format!("DELETE FROM edges WHERE id IN ({})", placeholders);

    // duckdb::params! requiere ownership, así que mapeamos los i64 a Value
    let params: Vec<&dyn rusqlite::ToSql> = edge_ids.iter().map(|id| id as &dyn rusqlite::ToSql).collect();

    // Usando execute con parámetros por posición
    conn.execute(sql.as_str(), edge_ids)?;
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

    let rows = stmt.query_map(params, |row| {
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

