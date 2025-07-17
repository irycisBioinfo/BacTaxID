use duckdb::{Connection, Result, params};
use serde::Deserialize;
use std::fs;

#[derive(Deserialize, Debug)]
pub struct MetadataConfig {
    pub genus: String,
    pub acronym: String,
    pub levels: String,
    pub kmer_size: i32,
    pub sketch_size: i32,
    pub click_size: i32,
    pub click_threshold: f64,
}

/// Estructura que encapsula la conexi�n a DuckDB.
pub struct DuckDb {
    conn: Connection,
}

impl DuckDb {
    /// Inicializa una nueva base de datos DuckDB en el archivo especificado.
    pub fn new(db_path: &str) -> Result<Self> {
        let conn = Connection::open(db_path)?;
        Ok(DuckDb { conn })
    }

    /// Crea la tabla `edges` con las columnas especificadas si no existe.
    pub fn init_edges_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS edges (
                id INTEGER PRIMARY KEY,
                source VARCHAR,
                target VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    /// Crea e inicializa la tabla `levels` con datos basados en el Vec<u64> proporcionado.
    /// Los valores de Level ser�n "L_0", "L_1", ..., "L_N" donde N es el n�mero de elementos.
    /// Los valores de Dist corresponden a los valores del Vec<u64> convertidos a DOUBLE.
    pub fn init_levels_table(&self, distances: Vec<u64>) -> Result<()> {
        // Crear la tabla si no existe
        let create_table_sql = "
            CREATE TABLE IF NOT EXISTS levels (
                level VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(create_table_sql)?;

        // Limpiar la tabla existente
        self.conn.execute("DELETE FROM levels", [])?;

        // Preparar la declaraci�n de inserci�n
        let mut stmt = self.conn.prepare("INSERT INTO levels (level, dist) VALUES (?, ?)")?;

        // Iterar sobre el vector e insertar cada elemento
        for (index, &distance) in distances.iter().enumerate() {
            let level_name = format!("L_{}", index);
            let dist_value = distance as f64;
            
            stmt.execute([&level_name, &dist_value.to_string()])?;
        }

        Ok(())
    }

    /// Crea din�micamente la tabla `code` basada en los niveles existentes en la tabla `levels`.
    /// La tabla tendr� una columna `sample` (VARCHAR, Primary Key) y columnas L_0 a L_N (INTEGER).
    pub fn init_code_table(&self) -> Result<()> {
        // Obtener los nombres de los niveles desde la tabla levels
        let mut stmt = self.conn.prepare("SELECT level FROM levels ORDER BY level")?;
        let mut rows = stmt.query([])?;
        
        let mut level_columns = Vec::new();
        while let Some(row) = rows.next()? {
            let level_name: String = row.get(0)?;
            level_columns.push(level_name);
        }

        // Construir la SQL para crear la tabla din�micamente
        let mut create_table_sql = String::from("CREATE TABLE IF NOT EXISTS code (\n    sample VARCHAR PRIMARY KEY");
        
        // Agregar columnas L_0 a L_N como INTEGER
        for level_name in &level_columns {
            create_table_sql.push_str(&format!(",\n    {} INTEGER", level_name));
        }
        
        create_table_sql.push_str("\n);");

        // Ejecutar la creaci�n de la tabla
        self.conn.execute_batch(&create_table_sql)?;
        
        Ok(())
    }


    /// Crea din�micamente la tabla `code_full` basada en los niveles existentes en la tabla `levels`.
    /// La tabla tendr� una columna `sample` (VARCHAR, Primary Key) y columnas L_0 a L_N (VARCHAR).
    pub fn init_code_full_table(&self) -> Result<()> {
        // Obtener los nombres de los niveles desde la tabla levels
        let mut stmt = self.conn.prepare("SELECT level FROM levels ORDER BY level")?;
        let mut rows = stmt.query([])?;
        
        let mut level_columns = Vec::new();
        while let Some(row) = rows.next()? {
            let level_name: String = row.get(0)?;
            level_columns.push(level_name);
        }

        // Construir la SQL para crear la tabla din�micamente
        let mut create_table_sql = String::from("CREATE TABLE IF NOT EXISTS code_full (\n    sample VARCHAR PRIMARY KEY");
        
        // Agregar columnas L_0 a L_N como VARCHAR
        for level_name in &level_columns {
            create_table_sql.push_str(&format!(",\n    {} VARCHAR", level_name));
        }
        
        create_table_sql.push_str("\n);");

        // Ejecutar la creaci�n de la tabla
        self.conn.execute_batch(&create_table_sql)?;
        
        Ok(())
    }

    /// Crea din�micamente la tabla `code_state` basada en los niveles existentes en la tabla `levels`.
    /// La tabla tendr� una columna `sample` (VARCHAR, Primary Key) y columnas L_0 a L_N (VARCHAR).
    pub fn init_code_state_table(&self) -> Result<()> {
        // Obtener los nombres de los niveles desde la tabla levels
        let mut stmt = self.conn.prepare("SELECT level FROM levels ORDER BY level")?;
        let mut rows = stmt.query([])?;
        
        let mut level_columns = Vec::new();
        while let Some(row) = rows.next()? {
            let level_name: String = row.get(0)?;
            level_columns.push(level_name);
        }

        // Construir la SQL para crear la tabla din�micamente
        let mut create_table_sql = String::from("CREATE TABLE IF NOT EXISTS code_state (\n    sample VARCHAR PRIMARY KEY");
        
        // Agregar columnas L_0 a L_N como VARCHAR
        for level_name in &level_columns {
            create_table_sql.push_str(&format!(",\n    {} VARCHAR", level_name));
        }
        
        create_table_sql.push_str("\n);");

        // Ejecutar la creaci�n de la tabla
        self.conn.execute_batch(&create_table_sql)?;
        
        Ok(())
    }

    /// Crea la tabla `metadata` con informaci�n sobre el experimento.
    /// La tabla contiene metadatos del an�lisis incluyendo genus, acronym, levels, kmer_size y sketch_size.
    pub fn init_metadata_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS metadata (
                genus VARCHAR,
                acronym VARCHAR,
                levels VARCHAR,
                kmer_size INTEGER,
                sketch_size INTEGER,
                click_size INTEGER,
                click_threshold DOUBLE
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }
    pub fn create_from_toml(toml_path: &str) -> Result<Self> {
        // 1. Leer y deserializar el archivo TOML primero para obtener el genus
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;

        // 2. Generar el nombre de la base de datos basado en el genus
        let db_name = format!("{}_bactaxid.db", config.genus);
        
        // 3. Crear una nueva instancia de DuckDb con el nombre generado
        let db = Self::new(&db_name)?;

        // 4. Inicializar todas las tablas usando los datos del TOML
        db.init_all_tables_from_config(&config)?;

        Ok(db)
    }

    /// Funci�n principal que inicializa toda la base de datos desde un archivo TOML.
    /// 
    /// # Argumentos
    /// 
    /// * `toml_path` - Ruta al archivo TOML que contiene la configuraci�n completa
    /// 
    /// # Ejemplo de archivo TOML
    /// 
    /// ```
    /// genus = "Escherichia"
    /// acronym = "EC"
    /// levels = "[0.95,0.98,0.99,0.999,0.9999]"
    /// kmer_size = 21
    /// sketch_size = 1000
    /// ```
    /// 
    /// # Funcionalidad
    /// 
    /// Esta funci�n realiza las siguientes acciones:
    /// 1. Lee y parsea el archivo TOML
    /// 2. Extrae y convierte el campo `levels` de string a Vec<f64>
    /// 3. Convierte los valores f64 a u64 multiplicando por 100000 (para preservar precisi�n)
    /// 4. Inicializa la tabla `metadata` con los datos del TOML
    /// 5. Inicializa la tabla `levels` con los valores parseados
    /// 6. Crea todas las dem�s tablas (edges, code, code_full, code_state)
    pub fn init_database_from_toml(&self, toml_path: &str) -> Result<()> {
        // 1. Leer y deserializar el archivo TOML
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;

        // 2. Parsear el campo levels usando funci�n asociada
        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e))?;

        // 3. Convertir Vec<f64> a Vec<u64>
        let levels_u64: Vec<u64> = levels_vec
            .iter()
            .map(|&x| (x * 100000.0) as u64)
            .collect();

        // 4. Inicializar tabla metadata
        self.init_metadata_table()?;
        self.insert_metadata_from_config(&config)?;

        // 5. Inicializar tabla levels
        self.init_levels_table(levels_u64)?;

        // 6. Inicializar todas las dem�s tablas
        self.init_edges_table()?;
        self.init_code_table()?;
        self.init_code_full_table()?;
        self.init_code_state_table()?;

        Ok(())
    }

    fn init_all_tables_from_config(&self, config: &MetadataConfig) -> Result<()> {
        // 2. Parsear el campo levels usando función asociada
        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e))?;

        // 3. Convertir Vec<f64> a Vec<u64>
        let levels_u64: Vec<u64> = levels_vec
            .iter()
            .map(|&x| (x * 100000.0) as u64)
            .collect();

        // 4. Inicializar tabla metadata
        self.init_metadata_table()?;
        self.insert_metadata_from_config(config)?;

        // 5. Inicializar tabla levels
        self.init_levels_table(levels_u64)?;

        // 6. Inicializar todas las demás tablas
        self.init_edges_table()?;
        self.init_code_table()?;
        self.init_code_full_table()?;
        self.init_code_state_table()?;

        Ok(())
    }
    
    /// Funci�n auxiliar para parsear el string de levels a Vec<f64>
    /// 
    /// # Argumentos
    /// 
    /// * `levels_str` - String en formato "[0.95,0.98,0.99,0.999,0.9999]"
    /// 
    /// # Retorna
    /// 
    /// Result<Vec<f64>, Box<dyn std::error::Error + Send + Sync>> - Vector de valores f64 o error de parsing
    fn parse_levels_string(levels_str: &str) -> Result<Vec<f64>, Box<dyn std::error::Error + Send + Sync>> {
        // Remover corchetes y espacios
        let cleaned = levels_str
            .trim()
            .strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .ok_or("El string de levels debe estar entre corchetes []")?;

        // Dividir por comas y parsear cada valor
        let levels: Result<Vec<f64>, _> = cleaned
            .split(',')
            .map(|s| s.trim().parse::<f64>())
            .collect();

        levels.map_err(|e| format!("Error parseando levels: {}", e).into())
    }

    fn insert_metadata_from_config(&self, config: &MetadataConfig) -> Result<()> {
        let mut stmt = self.conn.prepare(
            "INSERT INTO metadata (
                genus, acronym, levels, kmer_size, sketch_size, click_size, click_threshold
            ) VALUES (?, ?, ?, ?, ?, ?, ?)"
        )?;

        stmt.execute(params![
            config.genus,
            config.acronym,
            config.levels,
            config.kmer_size,
            config.sketch_size,
            config.click_size,
            config.click_threshold,
        ])?;

        Ok(())
    }

    /// Lee un archivo TOML y inserta los datos en la tabla metadata.
    /// Esta funci�n es independiente y puede ser usada por separado.
    /// 
    /// # Argumentos
    /// 
    /// * `toml_path` - Ruta al archivo TOML que contiene los metadatos
    /// 
    /// # Ejemplo de archivo TOML
    /// 
    /// ```
    /// genus = "Escherichia"
    /// acronym = "EC"
    /// levels = "[0.95,0.98,0.99,0.999,0.9999]"
    /// kmer_size = 21
    /// sketch_size = 1000
    /// ```
    pub fn load_metadata_from_toml(&self, toml_path: &str) -> Result<()> {
        // Leer el archivo TOML
        let toml_content = fs::read_to_string(toml_path)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;

        // Deserializar el contenido TOML
        let config: MetadataConfig = toml::from_str(&toml_content)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.to_string().into()))?;

        // Insertar datos usando la funci�n auxiliar
        self.insert_metadata_from_config(&config)?;

        Ok(())
    }

    /// Acceso a la conexi�n interna, si se requiere para operaciones avanzadas.
    pub fn connection(&self) -> &Connection {
        &self.conn
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    #[test]
    fn test_metadata_table_initialization() {
        let db = DuckDb::new(":memory:").expect("No se pudo crear la base de datos");
        db.init_metadata_table().expect("No se pudo crear la tabla metadata");

        // Verificar columnas nuevas en la tabla metadata
        let columns_check = "SELECT column_name FROM information_schema.columns WHERE table_name = 'metadata' ORDER BY ordinal_position";
        let mut stmt = db.conn.prepare(columns_check).expect("Error preparando consulta");
        let mut rows = stmt.query([]).expect("Error ejecutando consulta");
        let mut columns = Vec::new();
        while let Some(row) = rows.next().expect("Error leyendo fila") {
            let column_name: String = row.get(0).expect("Error obteniendo column_name");
            columns.push(column_name);
        }
        assert_eq!(
            columns,
            vec![
                "genus",
                "acronym",
                "levels",
                "kmer_size",
                "sketch_size",
                "click_size",
                "click_threshold"
            ]
        );
    }

    #[test]
    fn test_insert_and_read_metadata_with_click_fields() {
        let db = DuckDb::new(":memory:").expect("No se pudo crear la base de datos");
        db.init_metadata_table().expect("No se pudo crear la tabla metadata");

        let meta = MetadataConfig {
            genus: "Testus".to_string(),
            acronym: "TS".to_string(),
            levels: "[0.1,0.2]".to_string(),
            kmer_size: 15,
            sketch_size: 200,
            click_size: 42,
            click_threshold: 0.006,
        };
        db.insert_metadata_from_config(&meta).expect("No se pudo insertar");

        let mut stmt = db.conn.prepare(
            "SELECT genus, acronym, levels, kmer_size, sketch_size, click_size, click_threshold FROM metadata"
        ).unwrap();
        let mut rows = stmt.query([]).unwrap();
        if let Some(row) = rows.next().unwrap() {
            let genus: String = row.get(0).unwrap();
            let acronym: String = row.get(1).unwrap();
            let levels: String = row.get(2).unwrap();
            let kmer_size: i32 = row.get(3).unwrap();
            let sketch_size: i32 = row.get(4).unwrap();
            let click_size: i32 = row.get(5).unwrap();
            let click_threshold: f64 = row.get(6).unwrap();

            assert_eq!(genus, "Testus");
            assert_eq!(acronym, "TS");
            assert_eq!(levels, "[0.1,0.2]");
            assert_eq!(kmer_size, 15);
            assert_eq!(sketch_size, 200);
            assert_eq!(click_size, 42);
            assert!((click_threshold - 0.006).abs() < 1e-10);
        } else {
            panic!("No hay filas en metadata");
        }
    }

    #[test]
    fn test_all_from_toml_with_click_fields() {
        let db = DuckDb::new(":memory:").expect("No se pudo crear la base de datos");
        let toml_content = r#"
genus = "Clostridium"
acronym = "CL"
levels = "[0.9,0.95,0.99]"
kmer_size = 29
sketch_size = 1500
click_size = 77
click_threshold = 0.8
"#;

        let temp_file = "tmp_metadata_full.toml";
        {
            let mut file = fs::File::create(temp_file).unwrap();
            file.write_all(toml_content.as_bytes()).unwrap();
        }

        db.init_database_from_toml(temp_file).expect("Inicialización desde TOML falló");

        // Verificar que los nuevos campos están bien insertados en metadata
        let mut stmt = db.conn.prepare(
            "SELECT genus, acronym, levels, kmer_size, sketch_size, click_size, click_threshold FROM metadata"
        ).unwrap();
        let mut rows = stmt.query([]).unwrap();
        if let Some(row) = rows.next().unwrap() {
            assert_eq!(row.get::<_, String>(0).unwrap(), "Clostridium");
            assert_eq!(row.get::<_, String>(1).unwrap(), "CL");
            assert_eq!(row.get::<_, String>(2).unwrap(), "[0.9,0.95,0.99]");
            assert_eq!(row.get::<_, i32>(3).unwrap(), 29);
            assert_eq!(row.get::<_, i32>(4).unwrap(), 1500);
            assert_eq!(row.get::<_, i32>(5).unwrap(), 77);
       
        } else {
            panic!("No hay fila en metadata tras importación TOML");
        }

        fs::remove_file(temp_file).unwrap();
    }
}
