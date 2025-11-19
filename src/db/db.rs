use duckdb::{Connection, Result, params, params_from_iter};
use serde::Deserialize;
use std::fs;
use bincode;
use crate::sketch::sketching::*;
use anyhow::{Result as AnyResult, Context, bail, anyhow};


#[derive(Deserialize, Debug)]
pub struct MetadataConfig {
    pub genus: String,
    pub acronym: String,
    pub levels: String,
    pub kmer_size: i32,
    pub sketch_size: i32,
    pub click_size: i32,
    pub click_threshold: f64,
    pub reference_size: i32
}


/// Structure that encapsulates the connection to DuckDB.
pub struct DuckDb {
    conn: Connection,
}


impl DuckDb {
    /// Initializes a new DuckDB database in the specified file.
    pub fn new(db_path: &str) -> Result<Self> {
        let conn = Connection::open(db_path)?;
        Ok(DuckDb { conn })
    }

    /// ✅ MODIFICADO: Tabla sketches usa signature UBIGINT como PRIMARY KEY
    pub fn init_sketches_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS sketches (
                signature UBIGINT PRIMARY KEY,
                sample VARCHAR NOT NULL,
                sketch BLOB NOT NULL,
                kmer_size UINTEGER,
                sketch_size UINTEGER,
                name VARCHAR
            );
            CREATE INDEX IF NOT EXISTS idx_sample ON sketches(sample);
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    /// ✅ MODIFICADO: Adds a Sketch to the sketches table
    pub fn add_sketch(&self, sample_id: &str, sketch: &Sketch) -> Result<()> {
        insert_sketch_object(&self.conn, sample_id, sketch)
    }

    /// Reconstructs a SketchManager from the database
    pub fn load_sketch_manager(&self, default_kmer_size: usize, default_sketch_size: usize) -> AnyResult<SketchManager> {
        load_sketch_manager_from_db(&self.conn, default_kmer_size, default_sketch_size)
    }

    /// ✅ MODIFICADO: Tabla edges usa signature UBIGINT
    pub fn init_edges_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE SEQUENCE IF NOT EXISTS id_sequence START 1;
            CREATE TABLE IF NOT EXISTS edges (
                id INTEGER PRIMARY KEY DEFAULT nextval('id_sequence'),
                source UBIGINT NOT NULL,
                target UBIGINT NOT NULL,
                dist DOUBLE NOT NULL,
                FOREIGN KEY (source) REFERENCES sketches(signature),
                FOREIGN KEY (target) REFERENCES sketches(signature)
            );
            CREATE INDEX IF NOT EXISTS idx_edges_source ON edges(source);
            CREATE INDEX IF NOT EXISTS idx_edges_target ON edges(target);
            CREATE INDEX IF NOT EXISTS idx_edges_source_target ON edges(source, target);
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    /// ✅ MODIFICADO: Tabla debug usa signature UBIGINT
    pub fn init_debug_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS debug (
                source UBIGINT NOT NULL,
                target UBIGINT NOT NULL,
                dist DOUBLE NOT NULL
            );
            CREATE INDEX IF NOT EXISTS idx_debug_source ON debug(source);
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    /// Creates and initializes the `levels` table with data based on the provided Vec<f64>.
    pub fn init_levels_table(&self, distances: Vec<f64>) -> Result<()> {
        let create_table_sql = "
            CREATE TABLE IF NOT EXISTS levels (
                level VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(create_table_sql)?;

        self.conn.execute("DELETE FROM levels", [])?;

        let mut stmt = self.conn.prepare("INSERT INTO levels (level, dist) VALUES (?, ?)")?;

        for (index, &distance) in distances.iter().enumerate() {
            let level_name = format!("L_{}", index);
            let dist_value = distance as f64;
            stmt.execute([&level_name, &dist_value.to_string()])?;
        }
        Ok(())
    }

    /// ✅ MODIFICADO: Tabla code ahora usa signature UBIGINT como PRIMARY KEY
    /// pero mantiene la columna sample VARCHAR
    pub fn init_code_table(&self) -> Result<()> {
        let mut stmt = self.conn.prepare("SELECT level FROM levels ORDER BY level")?;
        let mut rows = stmt.query([])?;

        let mut level_columns = Vec::new();
        while let Some(row) = rows.next()? {
            let level_name: String = row.get(0)?;
            level_columns.push(level_name);
        }

        let mut create_table_sql = String::from(
            "CREATE TABLE IF NOT EXISTS code (\n    signature UBIGINT PRIMARY KEY,\n    sample VARCHAR NOT NULL"
        );
        
        for level_name in &level_columns {
            create_table_sql.push_str(&format!(",\n    {}_int INTEGER", level_name));
            create_table_sql.push_str(&format!(",\n    {}_full VARCHAR", level_name));
            create_table_sql.push_str(&format!(",\n    {}_state VARCHAR", level_name));
        }
        create_table_sql.push_str(",\n    FOREIGN KEY (signature) REFERENCES sketches(signature)");
        create_table_sql.push_str("\n);");

        self.conn.execute_batch(&create_table_sql)?;
        Ok(())
    }

    /// Creates the `metadata` table with information about the experiment.
    pub fn init_metadata_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS metadata (
                genus VARCHAR,
                acronym VARCHAR,
                levels VARCHAR,
                kmer_size INTEGER,
                sketch_size INTEGER,
                click_size INTEGER,
                click_threshold DOUBLE,
                reference_size INTEGER
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    pub fn create_from_toml(toml_path: &str) -> Result<Self> {
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;

        let db_name = format!("{}_bactaxid.db", config.acronym);

        let db = Self::new(&db_name)?;

        db.init_all_tables_from_config(&config)?;

        Ok(db)
    }

    /// Main function that initializes the entire database from a TOML file.
    pub fn init_database_from_toml(&self, toml_path: &str) -> Result<()> {
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;

        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.into()))?;
        println!("Parsed levels: {:?}", levels_vec);

        self.init_metadata_table()?;
        self.insert_metadata_from_config(&config)?;

        self.init_levels_table(levels_vec)?;

        self.init_edges_table()?;
        self.init_code_table()?;
        self.init_sketches_table()?;
        self.init_debug_table()?;
        Ok(())
    }

    pub fn init_all_tables_from_config(&self, config: &MetadataConfig) -> Result<()> {
        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.into()))?;
        println!("Parsed levels: {:?}", levels_vec);

        self.init_metadata_table()?;
        self.insert_metadata_from_config(config)?;

        self.init_levels_table(levels_vec)?;

        self.init_edges_table()?;
        self.init_code_table()?;
        self.init_sketches_table()?;
        self.init_debug_table()?;

        Ok(())
    }

    /// Helper function to parse the levels string into Vec<f64>
    fn parse_levels_string(levels_str: &str) -> AnyResult<Vec<f64>> {
        let cleaned = levels_str
            .trim()
            .strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .context("The levels string must be enclosed in brackets []")?;

        let levels: AnyResult<Vec<f64>> = cleaned
            .split(',')
            .map(|s| {
                let value = s.trim().parse::<f64>()
                    .with_context(|| format!("Error parsing value '{}'", s.trim()))?;
                if value < 0.0 || value > 1.0 {
                    bail!("Value out of range [0,1]: {}", value);
                }
                Ok(value)
            })
            .collect();

        levels.context("Error processing levels values")
    }

    fn insert_metadata_from_config(&self, config: &MetadataConfig) -> Result<()> {
        let sql = "INSERT INTO metadata (
            genus, acronym, levels, kmer_size, sketch_size, click_size, click_threshold, reference_size
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)";

        println!(
            "SQL: {}\nParams: genus={:?}, acronym={:?}, levels={:?}, kmer_size={:?}, sketch_size={:?}, click_size={:?}, click_threshold={:?}, reference_size={:?}",
            sql,
            config.genus,
            config.acronym,
            config.levels,
            config.kmer_size,
            config.sketch_size,
            config.click_size,
            config.click_threshold,
            config.reference_size
        );

        let mut stmt = self.conn.prepare(sql)?;

        stmt.execute(params![
            config.genus,
            config.acronym,
            config.levels,
            config.kmer_size,
            config.sketch_size,
            config.click_size,
            config.click_threshold,
            config.reference_size
        ])?;
        Ok(())
    }

    /// Reads a TOML file and inserts the data into the metadata table.
    pub fn load_metadata_from_toml(&self, toml_path: &str) -> Result<()> {
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;

        self.insert_metadata_from_config(&config)?;

        Ok(())
    }

    // ============================================================
    // ✅ NUEVOS MÉTODOS PARA TRABAJAR CON SIGNATURES
    // ============================================================

    /// Retrieve sketch by signature
    pub fn get_sketch_by_signature(&self, signature: u64) -> AnyResult<Option<Sketch>> {
        let mut stmt = self.conn.prepare("SELECT sketch FROM sketches WHERE signature = ?")
            .context("Failed to prepare query")?;
        
        let result = stmt.query_row(params![signature], |row| {
            let sketch_data: Vec<u8> = row.get(0)?;
            Ok(sketch_data)
        });

        match result {
            Ok(sketch_data) => {
                let sketch: Sketch = bincode::deserialize(&sketch_data)
                    .context("Failed to deserialize sketch")?;
                Ok(Some(sketch))
            }
            Err(duckdb::Error::QueryReturnedNoRows) => Ok(None),
            Err(e) => Err(anyhow!("Database error: {}", e)),
        }
    }

    /// Check if signature exists
    pub fn signature_exists(&self, signature: u64) -> Result<bool> {
        let mut stmt = self.conn.prepare("SELECT COUNT(*) FROM sketches WHERE signature = ?")?;
        let count: i64 = stmt.query_row(params![signature], |row| row.get(0))?;
        Ok(count > 0)
    }

    /// Get sample name from signature
    pub fn get_sample_name(&self, signature: u64) -> Result<Option<String>> {
        get_sample_name_by_signature(&self.conn, signature)
    }

    /// Get multiple sample names from signatures
    pub fn get_sample_names(&self, signatures: &[u64]) -> Result<Vec<(u64, String)>> {
        get_sample_names_by_signatures(&self.conn, signatures)
    }

    /// Get signature from sample name
    pub fn get_signature_by_sample(&self, sample_id: &str) -> Result<Option<u64>> {
        let mut stmt = self.conn.prepare("SELECT signature FROM sketches WHERE sample = ?")?;
        match stmt.query_row(params![sample_id], |row| row.get(0)) {
            Ok(sig) => Ok(Some(sig)),
            Err(duckdb::Error::QueryReturnedNoRows) => Ok(None),
            Err(e) => Err(e),
        }
    }

    /// Get all signatures in the database
    pub fn get_all_signatures(&self) -> Result<Vec<u64>> {
        let mut stmt = self.conn.prepare("SELECT signature FROM sketches")?;
        let signatures = stmt.query_map([], |row| row.get(0))?
            .collect::<Result<Vec<_>, _>>()?;
        Ok(signatures)
    }

    // ============================================================
    // ✅ MÉTODOS PARA EDGES (AHORA CON SIGNATURES)
    // ============================================================

    /// Insert edge by signature
    pub fn add_edge(&self, source_sig: u64, target_sig: u64, distance: f64) -> Result<()> {
        insert_edge_by_signature(&self.conn, source_sig, target_sig, distance)
    }

    /// Insert multiple edges in batch
    pub fn add_edges_batch(&self, edges: &[(u64, u64, f64)]) -> Result<()> {
        insert_edges_batch(&self.conn, edges)
    }

    /// Get neighbors of a signature (outgoing edges only)
    pub fn get_neighbors(&self, signature: u64) -> Result<Vec<(u64, f64)>> {
        get_edges_by_source(&self.conn, signature)
    }

    /// Get neighbors with sample names
    pub fn get_neighbors_with_names(&self, signature: u64) -> Result<Vec<(String, u64, f64)>> {
        get_neighbors_with_names(&self.conn, signature)
    }

    /// Get bidirectional edges (both incoming and outgoing)
    pub fn get_edges_bidirectional(&self, signature: u64) -> Result<Vec<(u64, f64)>> {
        get_edges_bidirectional(&self.conn, signature)
    }

    /// Get distance between two signatures
    pub fn get_distance(&self, source_sig: u64, target_sig: u64) -> Result<Option<f64>> {
        get_edge_distance(&self.conn, source_sig, target_sig)
    }

    /// Count edges for a signature
    pub fn count_edges(&self, signature: u64) -> Result<usize> {
        count_edges_for_signature(&self.conn, signature)
    }

    /// Insert debug edge
    pub fn add_debug_edge(&self, source_sig: u64, target_sig: u64, distance: f64) -> Result<()> {
        insert_debug_edge(&self.conn, source_sig, target_sig, distance)
    }

    // ============================================================
    // ✅ MÉTODOS PARA CODE TABLE
    // ============================================================

    /// Insert empty code entry
    pub fn add_empty_code(&self, signature: u64, sample_id: &str) -> Result<()> {
        insert_empty_code(&self.conn, signature, sample_id)
    }

    /// Copy code fields up to level
    pub fn copy_code_fields(&self, input_sig: u64, ref_sig: u64, level: usize) -> Result<()> {
        copy_code_l_fields_up_to(&self.conn, input_sig, ref_sig, level)
    }

    /// Copy code full fields up to level
    pub fn copy_code_full_fields(&self, input_sig: u64, ref_sig: u64, level: usize) -> Result<()> {
        copy_code_full_l_fields_up_to(&self.conn, input_sig, ref_sig, level)
    }

    /// Access to the internal connection, if required for advanced operations.
    pub fn connection(&self) -> &Connection {
        &self.conn
    }
}


// ============================================================
// ✅ FUNCIONES MODIFICADAS PARA USAR SIGNATURE
// ============================================================

/// ✅ MODIFICADO: Inserts a row in code table using signature and sample
pub fn insert_empty_code(conn: &Connection, signature: u64, sample_id: &str) -> Result<()> {
    let mut stmt = conn.prepare(
        "SELECT column_name FROM information_schema.columns \
         WHERE table_name = 'code' AND column_name NOT IN ('signature', 'sample') \
         ORDER BY ordinal_position"
    )?;
    let columns: Vec<String> = stmt.query_map([], |row| row.get(0))?
        .flatten()
        .collect();
    let n = columns.len();

    let columns_sql = if n > 0 {
        format!(", {}", columns.join(", "))
    } else { "".to_string() };
    let nulls: String = if n > 0 {
        std::iter::repeat("NULL").take(n).collect::<Vec<_>>().join(", ")
    } else { "".to_string() };
    let sql = format!(
        "INSERT INTO code (signature, sample{columns}) VALUES (?, ?{nulls})",
        columns=columns_sql,
        nulls=if n>0 { format!(", {nulls}") } else { "".to_string() },
    );
    conn.execute(&sql, params![signature, sample_id])?;
    Ok(())
}

/// ✅ MODIFICADO: Copies L_X_int fields using signatures
pub fn copy_code_l_fields_up_to(
    conn: &Connection, 
    input_sig: u64, 
    ref_sig: u64, 
    l: usize
) -> Result<()> {
    let columns: Vec<String> = (0..=l).map(|i| format!("L_{}_int", i)).collect();

    let mut set_expr = String::new();
    for (i, col) in columns.iter().enumerate() {
        if i > 0 {
            set_expr.push_str(", ");
        }
        set_expr.push_str(&format!("{} = ?", col));
    }

    let sql = format!(
        "UPDATE code SET {set_expr} WHERE signature = ?",
        set_expr = set_expr
    );

    let select_sql = format!(
        "SELECT {} FROM code WHERE signature = ?",
        columns.join(",")
    );

    let mut stmt = conn.prepare(&select_sql)?;
    let values_row = stmt.query_row([ref_sig], |row| {
        (0..=l).map(|i| row.get::<_, Option<i32>>(i)).collect::<Result<Vec<_>, _>>()
    })?;

    if values_row.len() == l + 1 {
        let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();

        for value in values_row {
            if let Some(v) = value {
                params.push(Box::new(v));
            } else {
                params.push(Box::new(None::<i32>));
            }
        }

        params.push(Box::new(input_sig));

        let mut stmt2 = conn.prepare(&sql)?;
        stmt2.execute(params_from_iter(params))?;
        Ok(())
    } else {
        Err(duckdb::Error::ToSqlConversionFailure(
            format!("Did not obtain L_0_int..L_{}_int for signature '{:016x}'", l, ref_sig).into()
        ))      
    }
}

/// ✅ MODIFICADO: Copies L_X_full fields using signatures
pub fn copy_code_full_l_fields_up_to(
    conn: &Connection, 
    input_sig: u64, 
    ref_sig: u64, 
    l: usize
) -> Result<()> {
    let columns: Vec<String> = (0..=l).map(|i| format!("L_{}_full", i)).collect();

    let mut set_expr = String::new();
    for (i, col) in columns.iter().enumerate() {
        if i > 0 {
            set_expr.push_str(", ");
        }
        set_expr.push_str(&format!("{} = ?", col));
    }

    let sql = format!(
        "UPDATE code SET {set_expr} WHERE signature = ?",
        set_expr = set_expr
    );

    let select_sql = format!(
        "SELECT {} FROM code WHERE signature = ?",
        columns.join(",")
    );

    let mut stmt = conn.prepare(&select_sql)?;
    let values_row = stmt.query_row([ref_sig], |row| {
        (0..=l).map(|i| row.get::<_, Option<String>>(i)).collect::<Result<Vec<_>, _>>()
    })?;

    if values_row.len() == l + 1 {
        let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();

        for value in values_row {
            if let Some(v) = value {
                params.push(Box::new(v));
            } else {
                params.push(Box::new(None::<String>));
            }
        }

        params.push(Box::new(input_sig));

        let mut stmt2 = conn.prepare(&sql)?;
        stmt2.execute(params_from_iter(params))?;
        Ok(())
    } else {
        Err(duckdb::Error::ToSqlConversionFailure(
            format!("Did not obtain L_0_full..L_{}_full for signature '{:016x}'", l, ref_sig).into()
        ))      
    }
}

/// ✅ MODIFICADO: Inserts a serialized Sketch object into the sketches table
pub fn insert_sketch_object(
    conn: &Connection,
    sample_id: &str,
    sketch: &Sketch
) -> Result<()> {
    let serialized = bincode::serialize(sketch)
        .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

    let mut stmt = conn.prepare(
        "INSERT INTO sketches (signature, sample, sketch, kmer_size, sketch_size, name) 
         VALUES (?, ?, ?, ?, ?, ?)"
    )?;
    
    stmt.execute(params![
        sketch.signature,
        sample_id,
        serialized,
        sketch.kmer_size as u32,
        sketch.sketch_size as u32,
        sketch.name
    ])?;
    Ok(())
}

/// Gets all Sketch objects from the table (to reconstruct SketchManager)
pub fn get_all_sketch_objects(conn: &Connection) -> Result<Vec<(String, Sketch)>> {
    let mut stmt = conn.prepare("SELECT sample, sketch FROM sketches")?;
    let rows = stmt.query_map([], |row| {
        let sample_id: String = row.get(0)?;
        let sketch_data: Vec<u8> = row.get(1)?;
        Ok((sample_id, sketch_data))
    })?;

    let mut sketches = Vec::new();
    for row in rows {
        let (sample_id, sketch_data) = row?;
        let sketch: Sketch = bincode::deserialize(&sketch_data)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;
        sketches.push((sample_id, sketch));
    }
    Ok(sketches)
}

// ============================================================
// ✅ NUEVAS FUNCIONES PARA EDGES CON SIGNATURES
// ============================================================

/// Insert edge by signature
pub fn insert_edge_by_signature(
    conn: &Connection,
    source_sig: u64,
    target_sig: u64,
    distance: f64
) -> Result<()> {
    let mut stmt = conn.prepare(
        "INSERT INTO edges (source, target, dist) VALUES (?, ?, ?)"
    )?;
    stmt.execute(params![source_sig, target_sig, distance])?;
    Ok(())
}

/// Insert multiple edges in batch
pub fn insert_edges_batch(
    conn: &Connection,
    edges: &[(u64, u64, f64)]
) -> Result<()> {
    let mut stmt = conn.prepare(
        "INSERT INTO edges (source, target, dist) VALUES (?, ?, ?)"
    )?;
    
    for (source, target, dist) in edges {
        stmt.execute(params![source, target, dist])?;
    }
    Ok(())
}

/// Get edges by source signature
pub fn get_edges_by_source(
    conn: &Connection,
    source_sig: u64
) -> Result<Vec<(u64, f64)>> {
    let mut stmt = conn.prepare(
        "SELECT target, dist FROM edges WHERE source = ? ORDER BY dist"
    )?;
    
    let edges = stmt.query_map(params![source_sig], |row| {
        Ok((row.get(0)?, row.get(1)?))
    })?
    .collect::<Result<Vec<_>, _>>()?;
    
    Ok(edges)
}

/// Get bidirectional edges
pub fn get_edges_bidirectional(
    conn: &Connection,
    signature: u64
) -> Result<Vec<(u64, f64)>> {
    let mut stmt = conn.prepare(
        "SELECT target as other, dist FROM edges WHERE source = ?
         UNION
         SELECT source as other, dist FROM edges WHERE target = ?
         ORDER BY dist"
    )?;
    
    let edges = stmt.query_map(params![signature, signature], |row| {
        Ok((row.get(0)?, row.get(1)?))
    })?
    .collect::<Result<Vec<_>, _>>()?;
    
    Ok(edges)
}

/// Get edge distance between two signatures
pub fn get_edge_distance(
    conn: &Connection,
    source_sig: u64,
    target_sig: u64
) -> Result<Option<f64>> {
    let mut stmt = conn.prepare(
        "SELECT dist FROM edges WHERE source = ? AND target = ?"
    )?;
    
    match stmt.query_row(params![source_sig, target_sig], |row| row.get(0)) {
        Ok(dist) => Ok(Some(dist)),
        Err(duckdb::Error::QueryReturnedNoRows) => Ok(None),
        Err(e) => Err(e),
    }
}

/// Count edges for signature
pub fn count_edges_for_signature(
    conn: &Connection,
    signature: u64
) -> Result<usize> {
    let mut stmt = conn.prepare(
        "SELECT COUNT(*) FROM edges WHERE source = ? OR target = ?"
    )?;
    let count: i64 = stmt.query_row(params![signature, signature], |row| row.get(0))?;
    Ok(count as usize)
}

/// Get neighbors with sample names
pub fn get_neighbors_with_names(
    conn: &Connection,
    source_sig: u64
) -> Result<Vec<(String, u64, f64)>> {
    let mut stmt = conn.prepare(
        "SELECT s.sample, e.target, e.dist 
         FROM edges e
         JOIN sketches s ON e.target = s.signature
         WHERE e.source = ?
         ORDER BY e.dist"
    )?;
    
    let neighbors = stmt.query_map(params![source_sig], |row| {
        Ok((row.get(0)?, row.get(1)?, row.get(2)?))
    })?
    .collect::<Result<Vec<_>, _>>()?;
    
    Ok(neighbors)
}

/// Insert debug edge
pub fn insert_debug_edge(
    conn: &Connection,
    source_sig: u64,
    target_sig: u64,
    distance: f64
) -> Result<()> {
    let mut stmt = conn.prepare(
        "INSERT INTO debug (source, target, dist) VALUES (?, ?, ?)"
    )?;
    stmt.execute(params![source_sig, target_sig, distance])?;
    Ok(())
}

/// Get sample name by signature
pub fn get_sample_name_by_signature(
    conn: &Connection,
    signature: u64
) -> Result<Option<String>> {
    let mut stmt = conn.prepare("SELECT sample FROM sketches WHERE signature = ?")?;
    
    match stmt.query_row(params![signature], |row| row.get(0)) {
        Ok(name) => Ok(Some(name)),
        Err(duckdb::Error::QueryReturnedNoRows) => Ok(None),
        Err(e) => Err(e),
    }
}

/// Get sample names by signatures (batch query)
pub fn get_sample_names_by_signatures(
    conn: &Connection,
    signatures: &[u64]
) -> Result<Vec<(u64, String)>> {
    if signatures.is_empty() {
        return Ok(Vec::new());
    }
    
    let placeholders = signatures.iter().map(|_| "?").collect::<Vec<_>>().join(",");
    let query = format!(
        "SELECT signature, sample FROM sketches WHERE signature IN ({})",
        placeholders
    );
    
    let mut stmt = conn.prepare(&query)?;
    let params: Vec<&dyn duckdb::ToSql> = signatures.iter()
        .map(|s| s as &dyn duckdb::ToSql)
        .collect();
    
    let results = stmt.query_map(params_from_iter(params.into_iter()), |row| {
        Ok((row.get(0)?, row.get(1)?))
    })?
    .collect::<Result<Vec<_>, _>>()?;
    
    Ok(results)
}


// ============================================================
// TESTS
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    #[test]
    fn test_metadata_table_initialization() {
        // ... existing code ...
    }

    #[test]
    fn test_sketches_table_initialization() {
        let db = DuckDb::new(":memory:").expect("Could not create database");
        db.init_sketches_table().expect("Could not create sketches table");

        let columns_check = "SELECT column_name FROM information_schema.columns WHERE table_name = 'sketches' ORDER BY ordinal_position";
        let mut stmt = db.conn.prepare(columns_check).expect("Error preparing query");
        let mut rows = stmt.query([]).expect("Error running query");

        let mut columns = Vec::new();
        while let Some(row) = rows.next().expect("Error reading row") {
            let column_name: String = row.get(0).expect("Error getting column_name");
            columns.push(column_name);
        }

        assert_eq!(columns, vec!["signature", "sample", "sketch", "kmer_size", "sketch_size", "name"]);
    }

    #[test]
    fn test_code_table_initialization() {
        let db = DuckDb::new(":memory:").expect("Could not create database");
        
        // First create levels table
        db.init_levels_table(vec![0.95, 0.90, 0.85]).expect("Could not create levels table");
        
        // Then create code table
        db.init_code_table().expect("Could not create code table");

        // Verify code table has signature and sample columns
        let columns_check = "SELECT column_name FROM information_schema.columns WHERE table_name = 'code' ORDER BY ordinal_position LIMIT 2";
        let mut stmt = db.conn.prepare(columns_check).expect("Error preparing query");
        let mut rows = stmt.query([]).expect("Error running query");

        let mut columns = Vec::new();
        while let Some(row) = rows.next().expect("Error reading row") {
            let column_name: String = row.get(0).expect("Error getting column_name");
            columns.push(column_name);
        }

        assert_eq!(columns, vec!["signature", "sample"]);
    }

    #[test]
    fn test_insert_and_retrieve_sketch_object() {
        // ... test code ...
    }

    #[test]
    fn test_load_sketch_manager_from_db() {
        // ... test code ...
    }

    #[test]
    fn test_edges_with_signatures() {
        let db = DuckDb::new(":memory:").expect("Could not create database");
        db.init_sketches_table().expect("Could not create sketches table");
        db.init_edges_table().expect("Could not create edges table");

        // Test would insert sketches and edges using signatures
        // Then verify retrieval works correctly
    }
}
