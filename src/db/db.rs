use duckdb::{Connection, Result, params, params_from_iter};
use serde::Deserialize;
use std::fs;
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


    /// Dynamically creates the sketches table based on sketch_size from metadata table
    pub fn init_sketches_table(&self) -> Result<()> {
        // Get sketch_size from metadata table
        let sketch_size: i32 = self.conn.query_row(
            "SELECT sketch_size FROM metadata LIMIT 1",
            [],
            |row| row.get(0)
        )?;

        // Build SQL with dynamic array size
        let schema_sql = format!(
            "CREATE TABLE IF NOT EXISTS sketches (
                signature UBIGINT PRIMARY KEY,
                name VARCHAR,
                sketch UBIGINT[{}]
            );",
            sketch_size
        );
        
        self.conn.execute_batch(&schema_sql)?;
        Ok(())
    }


    /// Adds a Sketch to the sketches table
    pub fn add_sketch(&self, sketch: &Sketch) -> Result<()> {
        insert_sketch_object(&self.conn, sketch)
    }


    /// Reconstructs a SketchManager from the database
    pub fn load_sketch_manager(&self, default_kmer_size: usize, default_sketch_size: usize) -> AnyResult<SketchManager> {
        load_sketch_manager_from_db(&self.conn, default_kmer_size, default_sketch_size)
    }


    /// Creates the `edges` table with the specified columns if it does not exist.
    pub fn init_edges_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE SEQUENCE IF NOT EXISTS id_sequence START 1;
            CREATE TABLE IF NOT EXISTS edges (
                id INTEGER PRIMARY KEY DEFAULT nextval('id_sequence'),
                source UBIGINT,
                target UBIGINT,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }


    /// Creates the `debug` table with the columns `Source`, `Target`, and `dist`.
    pub fn init_debug_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS debug (
                Source VARCHAR,
                Target VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }


    /// Creates the `duplicates` table with signature and sample columns.
    pub fn init_duplicates_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS duplicates (
                signature UBIGINT REFERENCES code(signature),
                sample VARCHAR
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }


    /// Creates and initializes the `levels` table with data based on the provided Vec<f64>.
    pub fn init_levels_table(&self, distances: Vec<f64>) -> Result<()> {
        // Create the table if it does not exist
        let create_table_sql = "
            CREATE TABLE IF NOT EXISTS levels (
                level VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(create_table_sql)?;


        // Clean the existing table
        self.conn.execute("DELETE FROM levels", [])?;


        // Prepare the insert statement
        let mut stmt = self.conn.prepare("INSERT INTO levels (level, dist) VALUES (?, ?)")?;


        // Iterate over the vector and insert each element
        for (index, &distance) in distances.iter().enumerate() {
            let level_name = format!("L_{}", index);
            let dist_value = distance as f64;
            stmt.execute([&level_name, &dist_value.to_string()])?;
        }
        Ok(())
    }


    /// Dynamically creates the merged `code` table based on the existing levels in the `levels` table.
    /// The table will have:
    /// - signature UBIGINT (Primary Key)
    /// - sample VARCHAR (sample name for queries and joins)
    /// - For each level L_X:
    ///   - L_X_int (INTEGER)
    ///   - L_X_full (VARCHAR)
    ///   - L_X_state (VARCHAR)
    pub fn init_code_table(&self) -> Result<()> {
        // Get the names of the levels from the levels table
        let mut stmt = self.conn.prepare("SELECT level FROM levels ORDER BY level")?;
        let mut rows = stmt.query([])?;


        let mut level_columns = Vec::new();
        while let Some(row) = rows.next()? {
            let level_name: String = row.get(0)?;
            level_columns.push(level_name);
        }


        // Build the SQL to create the table dynamically
        let mut create_table_sql = String::from(
            "CREATE TABLE IF NOT EXISTS code (\n    signature UBIGINT PRIMARY KEY,\n    sample VARCHAR"
        );
        
        // Add L_X_int, L_X_full, L_X_state columns for each level
        for level_name in &level_columns {
            create_table_sql.push_str(&format!(",\n    {}_int INTEGER", level_name));
            create_table_sql.push_str(&format!(",\n    {}_full VARCHAR", level_name));
            create_table_sql.push_str(&format!(",\n    {}_state VARCHAR", level_name));
        }
        create_table_sql.push_str("\n);");


        // Execute the table creation
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
        // 1. Read and deserialize the TOML file first to get the genus
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;


        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;


        // 2. Generate the database name based on the genus
        let db_name = format!("{}_bactaxid.db", config.acronym);


        // 3. Create a new instance of DuckDb with the generated name
        let db = Self::new(&db_name)?;


        // 4. Initialize all tables using the TOML data
        db.init_all_tables_from_config(&config)?;


        Ok(db)
    }


    /// Main function that initializes the entire database from a TOML file.
    pub fn init_database_from_toml(&self, toml_path: &str) -> Result<()> {
        // 1. Read and deserialize the TOML file
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;


        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;


        // 2. Parse the levels field using an associated function
        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.into()))?;
        println!("Parsed levels: {:?}", levels_vec);


        // 3. Initialize metadata table FIRST (required for sketches table)
        self.init_metadata_table()?;
        self.insert_metadata_from_config(&config)?;


        // 4. Initialize levels table
        self.init_levels_table(levels_vec)?;


        // 5. Initialize code table (depends on levels)
        self.init_code_table()?;
        
        // 6. Initialize other tables (sketches depends on metadata)
        self.init_edges_table()?;
        self.init_sketches_table()?;  // Now reads sketch_size from metadata
        self.init_duplicates_table()?;
        self.init_debug_table()?;
        
        Ok(())
    }


    pub fn init_all_tables_from_config(&self, config: &MetadataConfig) -> Result<()> {
        // 1. Parse the levels field using the associated function
        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.into()))?;
        println!("Parsed levels: {:?}", levels_vec);


        // 2. Initialize metadata table FIRST (required for sketches table)
        self.init_metadata_table()?;
        self.insert_metadata_from_config(config)?;


        // 3. Initialize levels table
        self.init_levels_table(levels_vec)?;


        // 4. Initialize code table (depends on levels)
        self.init_code_table()?;
        
        // 5. Initialize other tables (sketches depends on metadata)
        self.init_edges_table()?;
        self.init_sketches_table()?;  // Now reads sketch_size from metadata
        self.init_duplicates_table()?;
        self.init_debug_table()?;


        Ok(())
    }


    /// Helper function to parse the levels string into Vec<f64>
    fn parse_levels_string(levels_str: &str) -> AnyResult<Vec<f64>> {
        // Remove brackets and spaces
        let cleaned = levels_str
            .trim()
            .strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .context("The levels string must be enclosed in brackets []")?;


        // Split by commas and parse each value
        let levels: AnyResult<Vec<f64>> = cleaned
            .split(',')
            .map(|s| {
                let value = s.trim().parse::<f64>()
                    .with_context(|| format!("Error parsing value '{}'", s.trim()))?;
                // Validation: ensure it is between 0 and 1
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


        // Print the query and parameters for debugging
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
        // Read the TOML file
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;


        // Deserialize the TOML content
        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;


        // Insert data using the helper function
        self.insert_metadata_from_config(&config)?;


        Ok(())
    }


    /// Access to the internal connection, if required for advanced operations.
    pub fn connection(&self) -> &Connection {
        &self.conn
    }
}


/// Inserts a row in the merged code table
/// with the signature field equal to `signature`, sample equal to `sample_name`,
/// and the rest of the columns as NULL.
pub fn insert_empty_code(conn: &Connection, signature: u64, sample_name: &str) -> Result<()> {
    // Get all columns except the primary (signature) and sample
    let mut stmt = conn.prepare(
        "SELECT column_name FROM information_schema.columns \
         WHERE table_name = 'code' AND column_name NOT IN ('signature', 'sample') \
         ORDER BY ordinal_position"
    )?;
    let columns: Vec<String> = stmt.query_map([], |row| row.get(0))?
        .flatten()
        .collect();
    let n = columns.len();


    // Build the SQL: INSERT INTO code (signature, sample, col1, col2, ...) VALUES (?, ?, NULL, NULL, ...)
    let columns_sql = if n > 0 {
        format!(", {}", columns.join(", "))
    } else { "".to_string() };
    
    let nulls: String = if n > 0 {
        std::iter::repeat("NULL").take(n).collect::<Vec<_>>().join(", ")
    } else { "".to_string() };
    
    let sql = format!(
        "INSERT INTO code (signature, sample{columns}) VALUES (?, ?{nulls})",
        columns = columns_sql,
        nulls = if n > 0 { format!(", {}", nulls) } else { "".to_string() },
    );
    
    conn.execute(&sql, params![signature as u64, sample_name])?;
    Ok(())
}


/// Copies L_X_int fields from ref_sig to input_sig up to level l (inclusive)
pub fn copy_code_l_fields_up_to(conn: &Connection, input_sig: u64, ref_sig: u64, l: usize) -> Result<()> {
    // Generate list of columns L_0_int to L_L_int
    let columns: Vec<String> = (0..=l).map(|i| format!("L_{}_int", i)).collect();


    // Build UPDATE SQL
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


    // Retrieve values from ref_sig
    let select_sql = format!(
        "SELECT {} FROM code WHERE signature = ?",
        columns.join(",")
    );


    let mut stmt = conn.prepare(&select_sql)?;
    let values_row = stmt.query_row([ref_sig as u64], |row| {
        (0..=l).map(|i| row.get::<_, Option<i32>>(i)).collect::<Result<Vec<_>, _>>()
    })?;


    // Perform UPDATE using params_from_iter
    if values_row.len() == l + 1 {
        // Create parameter vector
        let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();


        // Add all L_0_int..L_L_int values
        for value in values_row {
            if let Some(v) = value {
                params.push(Box::new(v));
            } else {
                params.push(Box::new(None::<i32>));
            }
        }


        // Add input_sig at the end
        params.push(Box::new(input_sig as u64));


        // Execute using params_from_iter
        let mut stmt2 = conn.prepare(&sql)?;
        stmt2.execute(params_from_iter(params))?;
        Ok(())
    } else {
        Err(duckdb::Error::ToSqlConversionFailure(
            format!("Did not obtain L_0_int..L_{}_int for signature '{}'", l, ref_sig).into()
        ))      
    }
}


/// Copies L_X_full fields from ref_sig to input_sig up to level l (inclusive)
pub fn copy_code_full_l_fields_up_to(conn: &Connection, input_sig: u64, ref_sig: u64, l: usize) -> Result<()> {
    // Generate list of columns L_0_full to L_L_full
    let columns: Vec<String> = (0..=l).map(|i| format!("L_{}_full", i)).collect();


    // Build UPDATE SQL
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


    // Retrieve values from ref_sig
    let select_sql = format!(
        "SELECT {} FROM code WHERE signature = ?",
        columns.join(",")
    );


    let mut stmt = conn.prepare(&select_sql)?;
    let values_row = stmt.query_row([ref_sig as u64], |row| {
        (0..=l).map(|i| row.get::<_, Option<String>>(i)).collect::<Result<Vec<_>, _>>()
    })?;


    // Perform UPDATE using params_from_iter
    if values_row.len() == l + 1 {
        // Create parameter vector
        let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();


        // Add all L_0_full..L_L_full values
        for value in values_row {
            if let Some(v) = value {
                params.push(Box::new(v));
            } else {
                params.push(Box::new(None::<String>));
            }
        }


        // Add input_sig at the end
        params.push(Box::new(input_sig as u64));


        // Execute using params_from_iter
        let mut stmt2 = conn.prepare(&sql)?;
        stmt2.execute(params_from_iter(params))?;
        Ok(())
    } else {
        Err(duckdb::Error::ToSqlConversionFailure(
            format!("Did not obtain L_0_full..L_{}_full for signature '{}'", l, ref_sig).into()
        ))      
    }
}




/// Inserts a Sketch into the sketches table using UBIGINT[] array
/// Inserts a Sketch into the sketches table using UBIGINT[] array
pub fn insert_sketch_object(
    conn: &Connection,
    sketch: &Sketch
) -> Result<()> {
    let sketch_array: Vec<u64> = sketch.hashes().to_vec();
    
    // Convertir array a bytes para almacenar como binario
    let array_bytes: Vec<u8> = sketch_array
        .iter()
        .flat_map(|&val| val.to_le_bytes().to_vec())
        .collect();
    
    // Guardar como BLOB y luego CAST a UBIGINT[]
    let mut stmt = conn.prepare(
        "INSERT INTO sketches (signature, name, sketch) VALUES (?, ?, ?)"
    )?;
    
    stmt.execute(duckdb::params![
        sketch.signature as u64,
        &sketch.name,
        array_bytes
    ])?;
    
    Ok(())
}



/* /// Gets all Sketch objects from the table (to reconstruct SketchManager)
pub fn get_all_sketch_objects(conn: &Connection) -> Result<Vec<(u64, Sketch)>> {
    let mut stmt = conn.prepare("SELECT signature, name, sketch FROM sketches")?;
    let rows = stmt.query_map([], |row| {
        let signature: u64 = row.get(0)?;
        let name: String = row.get(1)?;
        // UBIGINT[] se devuelve como Vec<u8> o similar, necesitamos convertirlo
        let sketch_data: Vec<u8> = row.get(2)?;
        Ok((signature, name, sketch_data))
    })?;

    let mut sketches = Vec::new();
    for row in rows {
        let (signature, name, sketch_data) = row?;
        
        // Convertir Vec<u8> a Vec<u64> (8 bytes por elemento)
        let sketch_array: Vec<u64> = sketch_data
            .chunks_exact(8)
            .map(|chunk| u64::from_le_bytes([
                chunk[0], chunk[1], chunk[2], chunk[3],
                chunk[4], chunk[5], chunk[6], chunk[7],
            ]))
            .collect();
        
        // Reconstruir Sketch object
        let sketch = Sketch {
            signature,
            name,
            sketch: sketch_array,
            kmer_size: 0,  // Estos se cargarán desde metadata si es necesario
            sketch_size: 0,
        };
        
        sketches.push((signature, sketch));
    }
    Ok(sketches)
} */



// At the end of your duckdb.rs file
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
        
        // First insert metadata (required for sketches table)
        db.init_metadata_table().expect("Could not create metadata table");
        db.conn.execute(
            "INSERT INTO metadata (sketch_size) VALUES (?)",
            params![4000]
        ).expect("Could not insert metadata");
        
        db.init_sketches_table().expect("Could not create sketches table");


        // Check that the sketches table was created correctly
        let columns_check = "SELECT column_name FROM information_schema.columns WHERE table_name = 'sketches' ORDER BY ordinal_position";
        let mut stmt = db.conn.prepare(columns_check).expect("Error preparing query");
        let mut rows = stmt.query([]).expect("Error running query");


        let mut columns = Vec::new();
        while let Some(row) = rows.next().expect("Error reading row") {
            let column_name: String = row.get(0).expect("Error getting column_name");
            columns.push(column_name);
        }


        assert_eq!(columns, vec!["signature", "name", "sketch"]);
    }


    #[test]
    fn test_insert_and_retrieve_sketch_object() {
        // ... test code ...
    }


    #[test]
    fn test_load_sketch_manager_from_db() {
        // ... test code ...
    }
}
