use clap::Args;
use anyhow::{Result, Context, anyhow, bail};
use hashbrown::HashMap;
use std::fs;
use std::path::Path;
use rayon::prelude::*;
use duckdb::Connection;

use crate::{

    sketch::sketching::*,
    db::db::DuckDb
};

/// Argumentos del subcomando update
#[derive(Args, Debug)]
pub struct UpdateArgs {
    /// Ruta a la base de datos DuckDB
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,

    /// Ruta al SketchManager serializado
    #[arg(long, required = true, value_name = "SKETCH_PATH")]
    pub sketch: String,

    /// Número de CPUs para paralelizar con rayon
    #[arg(long, value_name = "N_CPUS", default_value_t = 1)]
    pub cpus: usize,

    /// Archivo de texto plano con paths a archivos FASTA (uno por línea)
    #[arg(long, required = true, value_name = "FILES_LIST")]
    pub files: String,
}

/// Contexto con la información precargada que necesita todo el flujo `update`
pub struct UpdateCtx<'a> {
    /// Conexión mutable a DuckDB (re-utilizable en todas las fases)
    pub conn: &'a mut Connection,
    /// Tabla levels (L_0 … L_N → dist)
    pub levels_map: HashMap<String, f64>,
    /// Toda la fila de metadata en forma (columna → valor en UTF-8)
    pub metadata_map: HashMap<String, String>,
}

impl<'a> UpdateCtx<'a> {
    /// Accesos de conveniencia para metadata
    pub fn kmer_size(&self) -> usize {
        self.metadata_map["kmer_size"].parse().unwrap_or(21)
    }

    pub fn sketch_size(&self) -> usize {
        self.metadata_map["sketch_size"].parse().unwrap_or(1000)
    }

    pub fn click_size(&self) -> usize {
        self.metadata_map["click_size"].parse().unwrap_or(50)
    }

    pub fn click_threshold(&self) -> f64 {
        self.metadata_map["click_threshold"].parse().unwrap_or(0.025)
    }

    pub fn genus(&self) -> &str {
        &self.metadata_map["genus"]
    }

    pub fn acronym(&self) -> &str {
        &self.metadata_map["acronym"]
    }

    pub fn levels(&self) -> &str {
        &self.metadata_map["levels"]
    }

    /// Acceso a la conexión mutable
    pub fn connection_mut(&mut self) -> &mut Connection {
        self.conn
    }
}

/// Lee la tabla `levels` y retorna un HashMap con los nombres de los niveles como clave
/// y la distancia (dist, f64) como valor.
fn load_levels_map(conn: &Connection) -> Result<HashMap<String, f64>> {
    let mut stmt = conn.prepare("SELECT level, dist FROM levels")
        .context("Error preparando consulta de tabla levels")?;
    
    let mut rows = stmt.query([])
        .context("Error ejecutando consulta de tabla levels")?;
    
    let mut map = HashMap::new();
    while let Some(row) = rows.next()? {
        let level: String = row.get(0)
            .context("Error obteniendo columna level")?;
        let dist: f64 = row.get(1)
            .context("Error obteniendo columna dist")?;
        map.insert(level, dist);
    }
    
    println!("✓ Cargados {} niveles desde la tabla levels", map.len());
    Ok(map)
}

/// Lee metadata (una sola fila) → HashMap<columna, valor_como_texto>
fn load_metadata_map(conn: &Connection) -> Result<HashMap<String, String>> {
    // Obtener las columnas dinámicamente
    let mut cols_stmt = conn.prepare(
        "SELECT column_name FROM information_schema.columns \
         WHERE table_name = 'metadata' ORDER BY ordinal_position"
    ).context("Error preparando consulta de columnas de metadata")?;
    
    let mut cols_rows = cols_stmt.query([])
        .context("Error ejecutando consulta de columnas")?;
    
    let mut cols = Vec::new();
    while let Some(row) = cols_rows.next()? {
        let col_name: String = row.get(0)
            .context("Error obteniendo nombre de columna")?;
        cols.push(col_name);
    }

    // Construir consulta SELECT dinámicamente
    let select_sql = format!("SELECT {} FROM metadata LIMIT 1", cols.join(", "));
    let mut row_stmt = conn.prepare(&select_sql)
        .context("Error preparando consulta de metadata")?;
    
    let row = row_stmt.query_row([], |r| {
        let mut map = HashMap::new();
        for (i, col) in cols.iter().enumerate() {
            // Todo metadata lo guardamos como texto para simplicidad
            let val: Option<String> = r.get(i)?;
            map.insert(col.clone(), val.unwrap_or_default());
        }
        Ok(map)
    }).context("Error leyendo fila de metadata")?;
    
    println!("✓ Cargados {} campos desde la tabla metadata", row.len());
    Ok(row)
}

/// Función principal del comando update
pub fn update_command(args: &UpdateArgs) -> Result<()> {
    println!("=== Iniciando comando update ===");
    println!("Database: {}", args.db);
    println!("SketchManager: {}", args.sketch);
    println!("CPUs: {}", args.cpus);
    println!("Lista de archivos: {}", args.files);

    // Verificar que los archivos existen
    if !Path::new(&args.db).exists() {
        return Err(anyhow::anyhow!("Base de datos no encontrada: {}", args.db));
    }
    if !Path::new(&args.sketch).exists() {
        return Err(anyhow::anyhow!("SketchManager no encontrado: {}", args.sketch));
    }
    if !Path::new(&args.files).exists() {
        return Err(anyhow::anyhow!("Lista de archivos no encontrada: {}", args.files));
    }

    // Abrir conexión mutable a la base de datos
    let mut conn = Connection::open(&args.db)
        .with_context(|| format!("Error abriendo base de datos: {}", args.db))?;

    println!("✓ Conexión a base de datos establecida");

    // Precargar tablas levels y metadata
    let levels_map = load_levels_map(&conn)?;
    let metadata_map = load_metadata_map(&conn)?;

    // Crear contexto compartido
    let mut ctx = UpdateCtx {
        conn: &mut conn,
        levels_map,
        metadata_map,
    };

    println!("\n=== Información del proyecto ===");
    println!("Genus: {}", ctx.genus());
    println!("Acronym: {}", ctx.acronym());
    println!("K-mer size: {}", ctx.kmer_size());
    println!("Sketch size: {}", ctx.sketch_size());
    println!("Click size: {}", ctx.click_size());
    println!("Click threshold: {}", ctx.click_threshold());
    println!("Niveles disponibles: {}", ctx.levels_map.len());

    // ✅ SOLUCIÓN: Extraer el valor del Result antes de usar
    let sketch_manager = SketchManager::load_from_disk(&args.sketch)
        .with_context(|| format!("Error cargando SketchManager: {}", args.sketch))?;

    println!("✓ SketchManager cargado. Contiene {} sketches", sketch_manager.length());



    // TODO: Aquí se agregará la lógica para:
    // 1. Comparar sketches query vs reference
    // 2. Actualizar contadores en level_counts
    // 3. Insertar/actualizar registros en tablas code, code_full, code_state
    // 4. Guardar SketchManager actualizado

    println!("\n=== Comando update completado exitosamente ===");
    Ok(())
}


/// Procesa **un** archivo FASTA:
/// 1. Crea un `Sketch` de la muestra usando `Sketch::new` (que lee el FASTA internamente).
/// 2. Crea el vector `query_class` con la longitud de los niveles.
/// 3. Devuelve `(sample_name, query_sketch, query_class)`
///
/// * `query_sketch` es el sketch creado desde el archivo FASTA.
/// * `query_class` es un vector de `u64` con longitud `ctx.levels_map.len()`
///   inicializado a 0; se usará posteriormente para contar coincidencias por nivel.
pub fn update_single_file(
    fasta_path: &Path,
    ctx: &mut UpdateCtx,
    _sketch_manager: &SketchManager, // No se modifica, solo se pasa por consistencia
) -> Result<()> {
    // -------- 1. Validaciones --------
    if !fasta_path.exists() {
        return Err(anyhow::anyhow!(
            "Archivo FASTA no encontrado: {:?}",
            fasta_path
        ));
    }

    // -------- 2. Crear el Sketch desde el path --------
    let kmer_size = ctx.kmer_size();
    let sketch_size = ctx.sketch_size();

    let mut query_sketch = SketchManager::new(ctx.kmer_size(), ctx.sketch_size());
    let new_sketch = Sketch::new(fasta_path,ctx.kmer_size(),ctx.sketch_size());

    query_sketch.add_sketch(new_sketch?)
    .with_context(|| format!("Error añadiendo sketch desde {:?}", fasta_path))?;

    // -------- 3. Extraer nombre de la muestra --------
    let sample_name = fasta_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    // -------- 4. Crear el vector de clasificación --------
    let levels_len = ctx.levels_map.len();
    let mut query_class: HashMap<_, usize> = ctx.levels_map
        .keys()  // Obtiene las claves
        .map(|key| (key.clone(), 0_usize))  // Mapea cada clave a (clave, 0)
        .collect();  // Recolecta en un HashMap



    println!(
        "✓ Sketch creado • muestra: {} • k: {} • sketch: {} • niveles: {}",
        sample_name,
        kmer_size,
        sketch_size,
        levels_len
    );

    let mut ref_db: Vec<String> = Vec::new();
    let mut distances: Vec<(String, String, f64)> = Vec::new();
    let mut best_hit: Option<(String, String, f64)> = None;

    let levels_data: Vec<_> = ctx.levels_map.iter()
    .map(|(k, v)| (k.clone(), v.clone()))
    .collect();

    for (level, Dist_level) in levels_data 
    { 
        if query_class.values().sum::<usize>() == 0 {
            ref_db = get_ids_by_level_value(&ctx.conn, "L_0", "C")?;
            
        } else {
            ref_db = find_matching_ids_by_levels_and_state(ctx.connection_mut(), &query_class, "C".to_string())?;
        }
        distances = pw_one_to_many(&query_sketch, &_sketch_manager, &ref_db);
        best_hit = find_best_hit(&distances);
        if let Some((sample, ref_name, dist)) = best_hit {
            if dist <= Dist_level {
                let value = get_level_value_by_id(ctx.connection_mut(), &ref_name, &level)?;
                query_class.insert(level.clone(), value.unwrap_or(0) as usize);
                // Aquí se actualizarían los contadores en query_class según corresponda
            } else {
                println!("No se encontró coincidencia cercana para {}", sample_name);
            }
        } else {
            println!("No se encontraron coincidencias para {}", sample_name);
        }
        
    }  

    Ok(())
}

pub fn get_level_value_by_id(
    conn: &Connection,
    id: &str,
    level_column: &str,
) -> Result<Option<usize>> {
    let column_name = format!("L_{}", level_column);
    let sql = format!("SELECT {} FROM code WHERE sample = ?", column_name);
    
    let mut stmt = conn.prepare(&sql)?;
    let mut rows = stmt.query([id])?;
    
    if let Some(row) = rows.next()? {
        let value: Option<i64> = row.get(0)?;
        Ok(value.map(|v| v as usize)) // Conversión directa
    } else {
        Ok(None)
    }
}

pub fn find_best_hit(
    data: &[(String, String, f64)],
) -> Option<(String, String, f64)> {
    data.iter()
        .filter(|t| !t.2.is_nan())                      // descartar NaN
        .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap())  // comparar f64
        .cloned()                                      // ← mueve (clona) toda la tupla
}


/// Devuelve los IDs (sample) de filas de `code` que coinciden con el vector de entrada
/// en las columnas L_0..L_N (solo valores != 0), y que 
/// en `code_state` tienen "C" en la columna L_{state_level}.
///
/// - `conn`: conexión DuckDB
/// - `levels`: vector con valores de L_0..L_N (los 0 son comodín, no filtran)
/// - `state_level`: índice del nivel para code_state (ej: 3 para L_3)

pub fn find_matching_ids_by_levels_and_state(
    conn: &Connection,
    levels: &HashMap<String, usize>,  // <- Cambiado a HashMap
    state_level: String,
) -> Result<Vec<String>> {
    // Construir condiciones sólo para los niveles != 0
    let mut conditions = Vec::new();
    let mut params: Vec<&dyn duckdb::ToSql> = Vec::new();
    
    for (column_id, val) in levels.iter() {
        if *val != 0 {
            let col = format!("L_{}", column_id);  // <- Usar column_id directamente
            conditions.push(format!("code.{} = ?", col));
            params.push(val);  // val es &usize, que implementa ToSql
        }
    }

    if conditions.is_empty() {
        bail!("El hashmap de niveles no tiene ningún valor distinto de cero; no se puede construir consulta.");
    }

    // Condición de estado (que code_state.L_x = 'C')
    let state_col = format!("L_{}", state_level);
    conditions.push(format!("cs.{} = 'C'", state_col));

    let where_clause = conditions.join(" AND ");
    let sql = format!(
        "SELECT code.sample \
         FROM code \
         JOIN code_state cs ON code.sample = cs.sample \
         WHERE {}",
        where_clause
    );

    let mut stmt = conn.prepare(&sql)?;
    let mut rows = stmt.query(params.as_slice())?;

    let mut result = Vec::new();
    while let Some(row) = rows.next()? {
        result.push(row.get(0)?);
    }
    Ok(result)
}

/// Devuelve los IDs (sample) de las filas de la tabla 'code'
/// que coinciden con el vector de entrada en las columnas L_0, L_1, ..., L_N,
/// pero solo se filtra por aquellas posiciones cuyo valor es distinto de 0.
/// 
/// - `conn`: conexión a DuckDB
/// - `levels`: hashmap de valores para L_0..L_N

/// Busca los samples cuya columna dinámica L_x es distinta de cero.
/// `levels_map` es un HashMap que mapea el sufijo de la columna (String) a su valor i64.
/// Busca los samples cuya columna dinámica L_x es distinta de cero.
pub fn find_matching_ids_by_levels(
    conn: &Connection,
    levels_map: &HashMap<String, i64>,
) -> Result<Vec<String>> {
    let mut conditions = Vec::new();
    let mut params: Vec<&dyn duckdb::ToSql> = Vec::new();

    // Iterar sobre (&String, &i64)
    for (column_id, val) in levels_map.iter() {
        if *val != 0 {
            let col = format!("L_{}", column_id);
            conditions.push(format!("{} = ?", col));
            params.push(val);  // val: &i64 vive en levels_map
        }
    }

    if conditions.is_empty() {
        bail!(
            "El HashMap de niveles no tiene ningún valor distinto de cero; no se puede construir consulta."
        );
    }

    let where_clause = conditions.join(" AND ");
    let sql = format!("SELECT sample FROM code WHERE {}", where_clause);

    let mut stmt = conn.prepare(&sql)?;
    let mut rows = stmt.query(params.as_slice())?;

    let mut result = Vec::new();
    while let Some(row) = rows.next()? {
        result.push(row.get(0)?);
    }
    Ok(result)
}


/// Obtiene todos los IDs de code_state donde una columna específica tiene un valor específico
/// 
/// # Parámetros
/// * `conn` - Conexión a DuckDB
/// * `level` - Nombre de la columna de nivel (ej: "L_0", "L_1", etc.)
/// * `value` - Valor a buscar (ej: "C", "G", etc.)
pub fn get_ids_by_level_value(
    conn: &Connection, 
    level: &str, 
    value: &str
) -> Result<Vec<String>> {
    // Validar que el nivel existe
    let column_check = "SELECT column_name FROM information_schema.columns 
                       WHERE table_name = 'code_state' AND column_name = ?";
    
    let mut check_stmt = conn.prepare(column_check)?;
    let mut check_rows = check_stmt.query([level])?;
    
    if check_rows.next()?.is_none() {
        bail!("La columna {} no existe en la tabla code_state", level);
        // o, equivalente:
        // return Err(anyhow!("La columna {} no existe en la tabla code_state", level));
    }


    // Construir consulta dinámica (cuidado con SQL injection)
    let sql = format!("SELECT sample FROM code_state WHERE {} = ?", level);
    
    let mut stmt = conn.prepare(&sql)?;
    let mut rows = stmt.query([value])?;
    
    let mut ids = Vec::new();
    while let Some(row) = rows.next()? {
        let sample_id: String = row.get(0)?;
        ids.push(sample_id);
    }
    
    println!("✓ Encontrados {} IDs con {} = '{}'", ids.len(), level, value);
    Ok(ids)
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_load_levels_map() {
        let conn = Connection::open(":memory:").unwrap();
        
        // Crear tabla levels de prueba
        conn.execute_batch("CREATE TABLE levels (level VARCHAR, dist DOUBLE)").unwrap();
        conn.execute_batch("INSERT INTO levels VALUES ('L_0', 95000.0), ('L_1', 98000.0), ('L_2', 99000.0)").unwrap();
        
        let result = load_levels_map(&conn).unwrap();
        
        assert_eq!(result.len(), 3);
        assert_eq!(result["L_0"], 95000.0);
        assert_eq!(result["L_1"], 98000.0);
        assert_eq!(result["L_2"], 99000.0);
    }

    #[test]
    fn test_load_metadata_map() {
        let conn = Connection::open(":memory:").unwrap();
        
        // Crear tabla metadata de prueba
        conn.execute_batch(
            "CREATE TABLE metadata (genus VARCHAR, acronym VARCHAR, kmer_size INTEGER, sketch_size INTEGER)"
        ).unwrap();
        conn.execute_batch(
            "INSERT INTO metadata VALUES ('Escherichia', 'EC', 21, 1000)"
        ).unwrap();
        
        let result = load_metadata_map(&conn).unwrap();
        
        assert_eq!(result.len(), 4);
        assert_eq!(result["genus"], "Escherichia");
        assert_eq!(result["acronym"], "EC");
        assert_eq!(result["kmer_size"], "21");
        assert_eq!(result["sketch_size"], "1000");
    }

    #[test]
    fn test_update_ctx_getters() {
        let mut conn = Connection::open(":memory:").unwrap();
        
        let mut metadata_map = HashMap::new();
        metadata_map.insert("kmer_size".to_string(), "25".to_string());
        metadata_map.insert("sketch_size".to_string(), "2000".to_string());
        metadata_map.insert("click_size".to_string(), "75".to_string());
        metadata_map.insert("click_threshold".to_string(), "0.05".to_string());
        
        let ctx = UpdateCtx {
            conn: &mut conn,
            levels_map: HashMap::new(),
            metadata_map,
        };
        
        assert_eq!(ctx.kmer_size(), 25);
        assert_eq!(ctx.sketch_size(), 2000);
        assert_eq!(ctx.click_size(), 75);
        assert!((ctx.click_threshold() - 0.05).abs() < f64::EPSILON);
    }
}
