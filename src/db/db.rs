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

    pub fn init_sketches_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS sketches (
                sample VARCHAR PRIMARY KEY,
                sketch BLOB
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }
    /// Añade un Sketch a la tabla sketches
    pub fn add_sketch(&self, sample_id: &str, sketch: &Sketch) -> Result<()> {
        insert_sketch_object(&self.conn, sample_id, sketch)
    }


    /// Reconstruye un SketchManager desde la base de datos
    pub fn load_sketch_manager(&self, default_kmer_size: usize, default_sketch_size: usize) -> AnyResult<SketchManager> {
        load_sketch_manager_from_db(&self.conn, default_kmer_size, default_sketch_size)
    }

    /// Crea la tabla `edges` con las columnas especificadas si no existe.
    pub fn init_edges_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE SEQUENCE id_sequence START 1;
            CREATE TABLE IF NOT EXISTS edges (
                id INTEGER PRIMARY KEY DEFAULT nextval('id_sequence'),
                source VARCHAR,
                target VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    // Tabla de debug para almacenar pares de nodos y distancias
    /// Crea la tabla `debug` con las columnas `Source`, `Target` y `dist`.
    /// Esta tabla se utiliza para almacenar pares de nodos y sus distancias.
    /// La columna `Source` es de tipo VARCHAR, `Target` es de tipo VARCHAR
    /// y `dist` es de tipo DOUBLE.
    ///
    pub fn init_debug_table(&self) -> Result<()> {
        let schema_sql = "
            CREATE TABLE IF NOT EXISTS debug (
                Source VARCHAR ,
                Target VARCHAR,
                dist DOUBLE
            );
        ";
        self.conn.execute_batch(schema_sql)?;
        Ok(())
    }

    /// Crea e inicializa la tabla `levels` con datos basados en el Vec<u64> proporcionado.
    /// Los valores de Level ser�n "L_0", "L_1", ..., "L_N" donde N es el n�mero de elementos.
    /// Los valores de Dist corresponden a los valores del Vec<f64> convertidos a DOUBLE.
    pub fn init_levels_table(&self, distances: Vec<f64>) -> Result<()> {
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
                click_threshold DOUBLE,
                reference_size INTEGER,
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
        let db_name = format!("{}_bactaxid.db", config.acronym);
        
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
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.into()))?;
        println!("Parsed levels: {:?}", levels_vec);
    

        // 4. Inicializar tabla metadata
        self.init_metadata_table()?;
        self.insert_metadata_from_config(&config)?;

        // 5. Inicializar tabla levels
        self.init_levels_table(levels_vec)?;

        // 6. Inicializar todas las dem�s tablas
        self.init_edges_table()?;
        self.init_code_table()?;
        self.init_code_full_table()?;
        self.init_code_state_table()?;
        self.init_sketches_table()?; 
        self.init_debug_table()?;
        Ok(())
    }

    pub fn init_all_tables_from_config(&self, config: &MetadataConfig) -> Result<()> {
        // 2. Parsear el campo levels usando función asociada
        let levels_vec = Self::parse_levels_string(&config.levels)
            .map_err(|e| duckdb::Error::ToSqlConversionFailure(e.into()))?;
        println!("Parsed levels: {:?}", levels_vec);

        // 4. Inicializar tabla metadata
        self.init_metadata_table()?;
        self.insert_metadata_from_config(config)?;

        // 5. Inicializar tabla levels
        self.init_levels_table(levels_vec)?;

        // 6. Inicializar todas las demás tablas
        self.init_edges_table()?;
        self.init_code_table()?;
        self.init_code_full_table()?;
        self.init_code_state_table()?;
        self.init_sketches_table()?;
        self.init_debug_table()?;

        Ok(())
    }
    


/// Función auxiliar para parsear el string de levels a Vec<f64>
/// 
/// # Argumentos
/// 
/// * `levels_str` - String en formato "[0.95,0.97,0.98,0.99,0.999]"
/// 
/// # Retorna
/// 
/// anyhow::Result<Vec<f64>> - Vector de valores f64 o error
fn parse_levels_string(levels_str: &str) -> AnyResult<Vec<f64>> {
    // Remover corchetes y espacios
    let cleaned = levels_str
        .trim()
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
        .context("El string de levels debe estar entre corchetes []")?;

    // Dividir por comas y parsear cada valor
    let levels: AnyResult<Vec<f64>> = cleaned
        .split(',')
        .map(|s| {
            let value = s.trim().parse::<f64>()
                .with_context(|| format!("Error parseando el valor '{}'", s.trim()))?;
            
            // Validación: asegurar que esté entre 0 y 1
            if value < 0.0 || value > 1.0 {
                bail!("Valor fuera del rango [0,1]: {}", value);
            }
            
            Ok(value)
        })
        .collect();

    levels.context("Error procesando los valores de levels")
}



 fn insert_metadata_from_config(&self, config: &MetadataConfig) -> Result<()> {
    let sql = "INSERT INTO metadata (
        genus, acronym, levels, kmer_size, sketch_size, click_size, click_threshold, reference_size
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)";

    // Imprime la consulta y los parámetros para debug
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
    /// genus = "Prueba"
    /// acronym = "PRB"
    /// levels = "[0.95,0.97,0.98,0.99,0.999]"
    /// kmer_size = 21
    /// sketch_size = 1200
    /// click_size = 8
    /// click_threshold = 0.8
    /// reference_size = 100
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


/// Inserta una fila en la tabla code-like (code, code_full o code_state)
/// con el campo sample igual a `sample_id` y el resto de columnas como NULL.
/// Asume que el nombre de la columna ID es `sample` y el resto de las columnas pueden variar y están después.
pub fn insert_empty_entry_code_table(
    conn: &Connection,
    table_name: &str,
    sample_id: &str
) -> Result<()> {
    // Obtiene todas las columnas menos la primaria (sample)
    let mut stmt = conn.prepare(&format!(
        "SELECT column_name FROM information_schema.columns \
         WHERE table_name = ? AND column_name != 'sample' ORDER BY ordinal_position"
    ))?;
    let columns: Vec<String> = stmt.query_map([table_name], |row| row.get(0))?
        .flatten()
        .collect();
    let n = columns.len();

    // Construye el SQL: INSERT INTO table (sample, col1, col2, ...) VALUES (?, NULL, NULL, ...)
    let columns_sql = if n > 0 {
        format!(", {}", columns.join(", "))
    } else { "".to_string() };
    let nulls: String = if n > 0 {
        std::iter::repeat("NULL").take(n).collect::<Vec<_>>().join(", ")
    } else { "".to_string() };
    let sql = format!(
        "INSERT INTO {table} (sample{columns}) VALUES (?{nulls})",
        table=table_name,
        columns=columns_sql,
        nulls=if n>0 { format!(", {nulls}") } else { "".to_string() },
    );
    // println!("SQL: {sql}");
    conn.execute(&sql, params![sample_id])?;
    Ok(())
}

pub fn insert_empty_code(conn: &Connection, sample_id: &str) -> Result<()> {
    insert_empty_entry_code_table(conn, "code", sample_id)
}
pub fn insert_empty_code_full(conn: &Connection, sample_id: &str) -> Result<()> {
    insert_empty_entry_code_table(conn, "code_full", sample_id)
}
pub fn insert_empty_code_state(conn: &Connection, sample_id: &str) -> Result<()> {
    insert_empty_entry_code_table(conn, "code_state", sample_id)
}


pub fn copy_l_fields_up_to(
    conn: &Connection,
    table_name: &str,
    input_id: &str,
    ref_id: &str,
    l: usize,
) -> Result<()> {
    // Genera lista de columnas L_0 a L_L
    let columns: Vec<String> = (0..=l).map(|i| format!("L_{}", i)).collect();
    
    // Construye SQL UPDATE
    let mut set_expr = String::new();
    for (i, col) in columns.iter().enumerate() {
        if i > 0 {
            set_expr.push_str(", ");
        }
        set_expr.push_str(&format!("{} = ?", col));
    }

    let sql = format!(
        "UPDATE {table} SET {set_expr} WHERE sample = ?",
        table = table_name,
        set_expr = set_expr
    );

    // Recupera los valores de ref_id
    let select_sql = format!(
        "SELECT {} FROM {} WHERE sample = ?",
        columns.join(","),
        table_name
    );
    
    let mut stmt = conn.prepare(&select_sql)?;
    let values_row = stmt.query_row([ref_id], |row| {
        (0..=l).map(|i| row.get::<_, Option<String>>(i)).collect::<Result<Vec<_>, _>>()
    })?;

    // Haz el UPDATE usando params_from_iter
    if values_row.len() == l + 1 {
        // Crear vector de parámetros
        let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();
        
        // Agregar todos los valores L_0..L_L
        for value in values_row {
            if let Some(v) = value {
                params.push(Box::new(v));
            } else {
                params.push(Box::new(None::<String>));
            }
        }
        
        // Agregar el input_id al final
        params.push(Box::new(input_id.to_string()));

        // Usar params_from_iter para ejecutar
        let mut stmt2 = conn.prepare(&sql)?;
        stmt2.execute(params_from_iter(params))?;
        
        Ok(())
    } else {
        Err(duckdb::Error::ToSqlConversionFailure(
            format!("No se obtuvieron L_0..L_{l} para sample '{ref_id}'").into()
        ))      
    }
}
pub fn copy_code_l_fields_up_to(conn: &Connection, input_id: &str, ref_id: &str, l: usize) -> Result<()> {
    copy_l_fields_up_to(conn, "code", input_id, ref_id, l)
}

pub fn copy_code_full_l_fields_up_to(conn: &Connection, input_id: &str, ref_id: &str, l: usize) -> Result<()> {
    copy_l_fields_up_to(conn, "code_full", input_id, ref_id, l)
}

/// Inserta un objeto Sketch serializado en la tabla sketches
pub fn insert_sketch_object(
    conn: &Connection,
    sample_id: &str,
    sketch: &Sketch
) -> Result<()> {
    let serialized = bincode::serialize(sketch)
        .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;
    
    let mut stmt = conn.prepare("INSERT INTO sketches (sample, sketch) VALUES (?, ?)")?;
    stmt.execute(params![sample_id, serialized])?;
    Ok(())
}


/// Obtiene todos los objetos Sketch de la tabla (para reconstruir SketchManager)
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

// Al final de tu archivo duckdb.rs
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;

    // Tus tests existentes...
    #[test]
    fn test_metadata_table_initialization() {
        // ... código existente
    }

    // Tests nuevos para sketches
    #[test]
    fn test_sketches_table_initialization() {
        let db = DuckDb::new(":memory:").expect("No se pudo crear la base de datos");
        db.init_sketches_table().expect("No se pudo crear la tabla sketches");

        // Verificar que la tabla sketches se creó correctamente
        let columns_check = "SELECT column_name FROM information_schema.columns WHERE table_name = 'sketches' ORDER BY ordinal_position";
        let mut stmt = db.conn.prepare(columns_check).expect("Error preparando consulta");
        let mut rows = stmt.query([]).expect("Error ejecutando consulta");
        
        let mut columns = Vec::new();
        while let Some(row) = rows.next().expect("Error leyendo fila") {
            let column_name: String = row.get(0).expect("Error obteniendo column_name");
            columns.push(column_name);
        }
        
        assert_eq!(columns, vec!["sample", "sketch"]);
    }

    #[test]
    fn test_insert_and_retrieve_sketch_object() {
        // ... test code
    }

    #[test]
    fn test_load_sketch_manager_from_db() {
        // ... test code
    }
}
