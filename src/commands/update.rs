use clap::Args;
use anyhow::{Result, Context, anyhow, bail};
use hashbrown::HashMap;
use std::fs;
use std::path::Path;
use std::time::Instant;
use rayon::prelude::*;
use duckdb::{Connection, Row, ToSql, params,Result as DuckResult};
use crate::graph::exists_clique;

use crate::{
    sketch::sketching::*,
    db::db::*,
    graph::*,
};

/// Argumentos del subcomando update
#[derive(Args, Debug)]
pub struct UpdateArgs {
    /// Ruta a la base de datos DuckDB
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,

    /// Número de CPUs para paralelizar con rayon
    #[arg(long, value_name = "N_CPUS", default_value_t = 1)]
    pub cpus: usize,

    /// Archivo de texto plano con paths a archivos FASTA (uno por línea)
    #[arg(long, required = true, value_name = "FILES_LIST")]
    pub files: String,

    #[arg(long, required= false, value_name = "DEBUG")]
    pub debug: bool
}

/// Contexto con la información precargada que necesita todo el flujo `update`
pub struct UpdateCtx<'a> {
    /// Conexión mutable a DuckDB (re-utilizable en todas las fases)
    pub conn: &'a mut Connection,
    /// Tabla levels (L_0 … L_N → dist) como vector ordenado
    pub levels: Vec<(String, f64)>,
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

    pub fn levels_str(&self) -> &str {
        &self.metadata_map["levels"]
    }

    /// Acceso a la conexión mutable
    pub fn connection_mut(&mut self) -> &mut Connection {
        self.conn
    }

    pub fn reference_size(&self) -> usize {
        self.metadata_map["reference_size"].parse().unwrap_or(100)
    }

    /// Número de niveles
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    /// Obtener distancia para un nivel específico
    pub fn get_level_distance(&self, level: &str) -> Option<f64> {
        self.levels.iter().find(|(l, _)| l == level).map(|(_, d)| *d)
    }
}

/// Lee la tabla `levels` y retorna un Vec<(String, f64)> ordenado por nivel
fn load_levels_vec(conn: &Connection) -> Result<Vec<(String, f64)>> {
    let start = Instant::now();
    
    let mut stmt = conn.prepare("SELECT level, dist FROM levels ORDER BY dist ASC")
        .context("Error preparando consulta de tabla levels")?;
    
    let mut rows = stmt.query([])
        .context("Error ejecutando consulta de tabla levels")?;
    
    let mut levels = Vec::new();
    while let Some(row) = rows.next()? {
        let level: String = row.get(0)
            .context("Error obteniendo columna level")?;
        let dist: f64 = row.get(1)
            .context("Error obteniendo columna dist")?;
        levels.push((level, dist));
    }
    
    let duration = start.elapsed();
    println!("✓ Cargados {} niveles desde la tabla levels [{:.2}ms]", 
             levels.len(), duration.as_millis());
    Ok(levels)
}

/// Lee metadata (una sola fila) → HashMap<columna, valor_como_texto>
fn load_metadata_map(conn: &Connection) -> Result<HashMap<String, String>> {
    let start = Instant::now();
    
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
            // Leer como Value genérico y convertir a String
            let val: Option<duckdb::types::Value> = r.get(i)?;
            let val_str = match val {
                Some(duckdb::types::Value::Int(i)) => i.to_string(),
                Some(duckdb::types::Value::Double(f)) => f.to_string(),
                Some(duckdb::types::Value::Text(s)) => s,
                Some(duckdb::types::Value::Boolean(b)) => b.to_string(),
                Some(v) => format!("{:?}", v), // Para otros tipos, usar debug format
                None => String::new(),
            };
            map.insert(col.clone(), val_str);
        }
        Ok(map)
    }).context("Error leyendo fila de metadata")?;
    
    let duration = start.elapsed();
    println!("✓ Cargados {} campos desde la tabla metadata [{:.2}ms]", 
             row.len(), duration.as_millis());
    Ok(row)
}

/// Función principal del comando update
pub fn update_command(args: &UpdateArgs) -> Result<()> {
    let command_start = Instant::now();
    
    println!("=== Iniciando comando update ===");
    println!("Database: {}", args.db);
    println!("CPUs: {}", args.cpus);
    println!("Lista de archivos: {}", args.files);

    // Verificar que los archivos existen
    let validation_start = Instant::now();
    if !Path::new(&args.db).exists() {
        return Err(anyhow::anyhow!("Base de datos no encontrada: {}", args.db));
    }
    if !Path::new(&args.files).exists() {
        return Err(anyhow::anyhow!("Lista de archivos no encontrada: {}", args.files));
    }
    println!("✓ Validación de archivos [{:.2}ms]", validation_start.elapsed().as_millis());

    // Abrir conexión mutable a la base de datos
    let db_connect_start = Instant::now();
    let mut conn = Connection::open(&args.db)
        .with_context(|| format!("Error abriendo base de datos: {}", args.db))?;
    println!("✓ Conexión a base de datos establecida [{:.2}ms]", 
             db_connect_start.elapsed().as_millis());

    // Precargar tablas levels y metadata
    let levels = load_levels_vec(&conn)?;
    let metadata_map = load_metadata_map(&conn)?;

    // Crear contexto compartido
    let mut ctx = UpdateCtx {
        conn: &mut conn,
        levels,
        metadata_map,
    };

    println!("\n=== Información del proyecto ===");
    println!("Genus: {}", ctx.genus());
    println!("Acronym: {}", ctx.acronym());
    println!("K-mer size: {}", ctx.kmer_size());
    println!("Sketch size: {}", ctx.sketch_size());
    println!("Click size: {}", ctx.click_size());
    println!("Click threshold: {}", ctx.click_threshold());
    println!("Niveles disponibles: {}", ctx.levels.len());

    // Cargar SketchManager desde DuckDB usando UpdateCtx
    let sketch_load_start = Instant::now();
    let mut sketch_manager_result = load_sketch_manager_from_db(
        ctx.conn,
        ctx.kmer_size(),
        ctx.sketch_size()
    );
    println!("✓ SketchManager carga desde DB [{:.2}ms]", 
             sketch_load_start.elapsed().as_millis());

    match &mut sketch_manager_result {
        Ok(sketch_manager) => {
            println!(
                "✓ SketchManager cargado desde DB. Contiene {} sketches",
                sketch_manager.length()
            );

            // Leer lista de archivos
            let file_read_start = Instant::now();
            let file_content = fs::read_to_string(&args.files)
                .with_context(|| format!("Error leyendo archivo de lista de archivos: {}", args.files))?;
            let files: Vec<&str> = file_content.lines().collect();
            println!("✓ Lista de archivos leída: {} archivos [{:.2}ms]", 
                     files.len(), file_read_start.elapsed().as_millis());

            // Procesar cada archivo
            let mut processed_files = 0;
            let mut skipped_files = 0;
            let processing_start = Instant::now();

            for line in files {
                let fasta_path = Path::new(line.trim());
                if !fasta_path.exists() {
                    eprintln!("Archivo FASTA no encontrado: {}", fasta_path.display());
                    continue;
                }

                let sample_name = fasta_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string();

                if sketch_manager.contains(&sample_name) {
                    println!("✓ El sketch para {} ya existe, se omitirá.", sample_name);
                    skipped_files += 1;
                    continue;
                }

                let file_start = Instant::now();
                update_single_file(fasta_path, &mut ctx, sketch_manager, args.debug)?;
                processed_files += 1;
                println!("  └─ Archivo procesado [{:.2}s]", file_start.elapsed().as_secs_f32());
            }

            println!("✓ Procesamiento de archivos completado: {} procesados, {} omitidos [{:.2}s]",
                     processed_files, skipped_files, processing_start.elapsed().as_secs_f32());
        }
        Err(ref e) => {
            eprintln!("Error al cargar SketchManager desde la base de datos: {e}");
        }
    }

    let total_duration = command_start.elapsed();
    println!("\n=== Comando update completado exitosamente [{:.2}s] ===", 
             total_duration.as_secs_f32());
    Ok(())
}

pub struct Query {
    pub sample_name: String,
    pub code: Vec<usize>,         // vector de enteros positivos
    pub code_full: Vec<String>,   // vector de strings
    pub code_state: Vec<String>,  // vector de strings
    pub sketch: SketchManager,
}

impl Query {
    /// Crea un nuevo Query cargando el SketchManager+Sketch, y vectores del tamaño adecuado
    pub fn new(path: &Path, ctx: &UpdateCtx) -> anyhow::Result<Self> {
        let query_creation_start = Instant::now();
        
        let sample_name = path.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
        let n_levels = ctx.levels.len();

        // Vector de enteros positivos (inicializados a 0)
        let code = vec![0_usize; n_levels];
        // code_full inicializado con strings vacíos ("")
        let code_full = vec!["".to_string(); n_levels];
        // code_state inicializado con "S" (sin clasificar por defecto)
        let code_state = vec!["S".to_string(); n_levels];
        
        // Instanciar SketchManager y cargar Sketch
        let sketch_creation_start = Instant::now();
        let mut sketch_manager = SketchManager::new(ctx.kmer_size(), ctx.sketch_size());
        let sketch = Sketch::new(path, ctx.kmer_size(), ctx.sketch_size())?;
        sketch_manager.add_sketch(sketch)?;
        println!("    ✓ Sketch creado [{:.2}ms]", sketch_creation_start.elapsed().as_millis());

        println!("    ✓ Query inicializado [{:.2}ms]", query_creation_start.elapsed().as_millis());

        Ok(Query {
            sample_name,
            code,
            code_full,
            code_state,
            sketch: sketch_manager,
        })
    }
}

/// Procesa **un** archivo FASTA
pub fn update_single_file(
    fasta_path: &Path,
    ctx: &mut UpdateCtx,
    _sketch_manager: &mut SketchManager,
    debug: bool
) -> Result<()> {
    let file_processing_start = Instant::now();
    
    // -------- 1. Crear Query (incluye validaciones, sketch y vectores) --------
    let query_start = Instant::now();
    let mut query = Query::new(fasta_path, ctx)?;
    println!("  ✓ Query creado • muestra: {} • k: {} • sketch: {} • niveles: {} [{:.2}ms]",
        query.sample_name,
        ctx.kmer_size(),
        ctx.sketch_size(),
        ctx.num_levels(),
        query_start.elapsed().as_millis()
    );

    // -------- 2. Variables para el procesamiento --------
    let mut ref_db: Vec<(String, usize, String, String)>;
    let mut distances: Vec<(String, String, f64)>;
    let mut bh: Option<(String, usize, String, String)> = None;
    let mut ref_ids: Vec<String> = Vec::new();

    let click_threshold = ctx.click_threshold();
    let click_size = ctx.click_size();
    let reference_size = ctx.reference_size();

    // -------- 3. Procesar cada nivel en orden ascendente de distancia --------
    for i in 0..ctx.levels.len() {
        let level_start = Instant::now();
        println!("    🔍 Procesando nivel {} [{:.4} threshold]", i, ctx.levels[i].1);
        
        // 3.1 Recuperar clasificadores iniciales (condition = "C")
        let classifiers_start = Instant::now();
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "C")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "C")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
        println!("      ✓ Clasificadores obtenidos: {} [{:.2}ms]", 
                 ref_ids.len(), classifiers_start.elapsed().as_millis());

        // 3.2 Calcular distancias
        let distances_start = Instant::now();
        distances = pw_one_to_many(
            &query.sketch,
            _sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        println!("      ✓ Distancias calculadas: {} [{:.2}ms]", 
                 distances.len(), distances_start.elapsed().as_millis());

        // --- DEBUG: guardar distances si debug == true ---
        if debug {
            let debug_start = Instant::now();
            let tx = ctx.connection_mut().transaction()?;
            for (q, r, d) in &distances {
                tx.execute(
                    "INSERT INTO debug (Source, Target,dist) VALUES (?1, ?2, ?3)",
                    params![q, r, d],
                )?;
            }
            tx.commit()?;
            println!("      ✓ Debug data guardado [{:.2}ms]", debug_start.elapsed().as_millis());
        }
        
        // 3.3 Procesar distancias
        if !distances.is_empty() {
            let best_hit_start = Instant::now();
            // 5. Actualizar según best_hit
            bh = best_hit(&distances, &ref_db);
            println!("      ✓ Best hit calculado [{:.2}ms]", best_hit_start.elapsed().as_millis());

            // Obtener valores de best_hit
            let code_val = bh.as_ref().map_or(0, |(_, code, _, _)| *code);
            let code_full_val = bh
                .as_ref()
                .map_or_else(|| "".to_string(), |(_, _, cf, _)| cf.clone());

            if code_val != 0 {
                // Actualizar query con best hit válido
                query.code[i] = code_val;
                query.code_full[i] = code_full_val.clone();

                if i > 0 && query.code_full[i].starts_with(query.code_full[i-1].as_str()) == false {
                    println!("Warning: El code_full del best hit '{}' no es consistente con el nivel anterior '{}' en nivel {}, tratando como no hay candidatos", 
                        query.code_full[i], query.code_full[i-1], i);
                    panic!("Inconsistencia en code_full entre niveles consecutivos");
                }

                let classifier_check_start = Instant::now();
                if is_classifier(
                    &distances,
                    bh.as_ref().unwrap(),
                    &ref_db,
                    click_threshold,
                    reference_size,
                ) {
                    query.code_state[i] = "C".to_string();
                }
                println!("      ✓ Clasificador evaluado [{:.2}ms]", classifier_check_start.elapsed().as_millis());
                
                println!("    ✓ Nivel {} completado con best hit [{:.2}ms]", i, level_start.elapsed().as_millis());
                // Continuar al siguiente nivel
                continue;
            } else {
                println!(
                    "Warning: Best hit code = 0 en nivel {}, tratando como no hay candidatos",
                    i
                );
            }
        }

        // 3.4 No hay candidatos válidos (o code == 0): buscar con condition = "ALL"
        println!(
            "      No hay candidatos válidos en nivel {}, buscando con condition = 'ALL'...",
            i
        );
        
        let all_classifiers_start = Instant::now();
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "ALL")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "ALL")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
        println!("      ✓ Referencias (ALL) obtenidas: {} [{:.2}ms]", 
                 ref_ids.len(), all_classifiers_start.elapsed().as_millis());

        let all_distances_start = Instant::now();
        distances = pw_one_to_many(
            &query.sketch,
            _sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        println!("      ✓ Distancias (ALL) calculadas: {} [{:.2}ms]", 
                 distances.len(), all_distances_start.elapsed().as_millis());

        // --- DEBUG: guardar distances si debug == true ---
        if debug {
            let debug_start = Instant::now();
            let tx = ctx.connection_mut().transaction()?;
            for (q, r, d) in &distances {
                tx.execute(
                    "INSERT INTO debug (Source, Target,dist) VALUES (?1, ?2, ?3)",
                    params![q, r, d],
                )?;
            }
            tx.commit()?;
            println!("      ✓ Debug data (ALL) guardado [{:.2}ms]", debug_start.elapsed().as_millis());
        }

        if !distances.is_empty() {
            // 3.5 Buscar cliques y crear nuevos grupos
            let cliques_start = Instant::now();
            let levels = ctx.levels.clone();
            let click_size = ctx.click_size();
            let conn = ctx.connection_mut();
            look_for_cliques(conn, &distances, i, &levels, click_size, &mut query)?;
            println!("      ✓ Búsqueda de cliques completada [{:.2}ms]", cliques_start.elapsed().as_millis());
            println!("    ✓ Nivel {} completado con cliques [{:.2}ms]", i, level_start.elapsed().as_millis());
            break;
        } else {
            println!("    ○ Nivel {} sin distancias válidas [{:.2}ms]", i, level_start.elapsed().as_millis());
            break;
        }
    }

    // -------- 4. Actualizar la base de datos y guardar sketch --------
    let db_update_start = Instant::now();
    update_duckdb(ctx.connection_mut(), &query)?;
    println!("    ✓ Base de datos actualizada [{:.2}ms]", db_update_start.elapsed().as_millis());

    let sketch_save_start = Instant::now();
    let sketch = query
        .sketch
        .get_sketch(&query.sample_name)
        .expect("El sketch debería existir después de Query::new()");
    insert_sketch_object(ctx.conn, &query.sample_name, sketch)
        .with_context(|| format!("Error guardando sketch para muestra: {}", query.sample_name))?;
    _sketch_manager
        .add_sketch(sketch.clone())
        .with_context(|| format!("Error añadiendo sketch para muestra: {}", query.sample_name))?;
    println!("    ✓ Sketch guardado [{:.2}ms]", sketch_save_start.elapsed().as_millis());

    println!(
        "✓ Procesamiento completo • muestra: {} • código final: {:?} [{:.2}s]",
        query.sample_name, query.code, file_processing_start.elapsed().as_secs_f32()
    );

    Ok(())
}

/// Recupera clasificadores de la tabla code fusionada
pub fn retrieve_classifiers(
    conn: &Connection,
    level: usize,
    group: &str,
    condition: &str,
) -> Result<Vec<(String, usize, String, String)>> {
    let query_start = Instant::now();
    
    let level_int_col = format!("L_{}_int", level);
    let level_full_col = format!("L_{}_full", level);
    let level_state_col = format!("L_{}_state", level);
    
    // Para los filtros, necesitamos usar level-1 (excepto cuando level = 0)
    let filter_full_col = if level == 0 {
        format!("L_{}_full", level)
    } else {
        format!("L_{}_full", level - 1)
    };
    
    let (sql, params): (String, Vec<&dyn duckdb::ToSql>) = match (level, condition) {
        // Level 0 con condición 'C'
        (0, "C") => (
            format!(
                "SELECT 
                    sample,
                    {level_int_col} as code,
                    {level_full_col} as code_full,
                    {level_state_col} as code_state
                FROM code 
                WHERE {level_state_col} = 'C'"
            ),
            Vec::new()
        ),
        
        // Level 0 con condición 'ALL'
        (0, "ALL") => (
            format!(
                "SELECT 
                    sample,
                    {level_int_col} as code,
                    {level_full_col} as code_full,
                    {level_state_col} as code_state
                FROM code 
                WHERE {level_full_col} = ''"
            ),
            Vec::new()
        ),
        
        // Level > 0 con condición 'C'
        (level, "C") if level > 0 => (
            format!(
                "SELECT 
                    sample,
                    {level_int_col} as code,
                    {level_full_col} as code_full,
                    {level_state_col} as code_state
                FROM code 
                WHERE {level_state_col} = 'C' 
                AND {filter_full_col} = ?"
            ),
            vec![&group as &dyn duckdb::ToSql]
        ),
        
        // Level > 0 con condición 'ALL'
        (level, "ALL") if level > 0 => (
            format!(
                "SELECT 
                    sample,
                    {level_int_col} as code,
                    {level_full_col} as code_full,
                    {level_state_col} as code_state
                FROM code 
                WHERE {filter_full_col} = ?
                AND {level_full_col} = ''"
            ),
            vec![&group as &dyn duckdb::ToSql]
        ),
        
        // Condición no reconocida
        _ => {
            bail!("Invalid condition '{}'. Use 'C' or 'ALL'", condition);
        }
    };

    let mut stmt = conn.prepare(&sql)?;
    let mut rows = stmt.query(params.as_slice())?;

    let mut result = Vec::new();
    let mut row_count = 0;
    
    while let Some(row) = rows.next()? {
        row_count += 1;
        
        // Manejo seguro de cada columna que puede ser NULL
        let sample: String = match row.get::<_, Option<String>>(0)? {
            Some(val) => val,
            None => {
                println!("WARNING: NULL sample at row {}, skipping", row_count);
                continue;
            }
        };
        
        let code: usize = match row.get::<_, Option<i64>>(1)? {
            Some(val) => {
                if val < 0 {
                    println!("WARNING: Negative code {} at row {}, skipping", val, row_count);
                    continue;
                }
                val as usize
            },
            None => {
                println!("WARNING: NULL code at row {}, skipping", row_count);
                continue;
            }
        };
        
        let code_full: String = match row.get::<_, Option<String>>(2)? {
            Some(val) => val,
            None => {
                println!("WARNING: NULL code_full at row {}, using 'UNKNOWN'", row_count);
                "UNKNOWN".to_string()
            }
        };
        
        let code_state: String = match row.get::<_, Option<String>>(3)? {
            Some(val) => val,
            None => {
                println!("WARNING: NULL code_state at row {}, using 'UNKNOWN'", row_count);
                "UNKNOWN".to_string()
            }
        };
        
        result.push((sample, code, code_full, code_state));
    }

    let duration = query_start.elapsed();
    if duration.as_millis() > 10 { // Solo reportar si toma más de 10ms
        println!("        └─ SQL query ejecutada [{:.2}ms]", duration.as_millis());
    }

    Ok(result)
}

/// Encuentra la similitud maxima y devuelve la tupla correspondiente de ref_db.
pub fn best_hit(
    distances: &[(String, String, f64)],
    ref_db: &[(String, usize, String, String)],
) -> Option<(String, usize, String, String)> {
    // 1. Encontrar la distancia mínima (ignorar NaN)
    let best_distance = distances
        .iter()
        .filter(|(_, _, dist)| !dist.is_nan())
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())?;
    
    // 2. Obtener el ref_name del mejor hit
    let best_ref_name = &best_distance.1;
    
    // 3. Buscar la tupla correspondiente en ref_db
    ref_db
        .iter()
        .find(|(sample, _, _, _)| sample == best_ref_name)
        .cloned()
}

/// Determina si un best_hit pasa el criterio de clasificación.
pub fn is_classifier(
    distances: &[(String, String, f64)],
    best_hit: &(String, usize, String, String),
    ref_db: &[(String, usize, String, String)],
    click_threshold: f64,
    reference_size: usize,
) -> bool {
    // Extraer el code_full del best_hit
    let best_full = &best_hit.2;

    // 1. Contar cuántas filas en ref_db tienen ese mismo code_full
    let count_ref = ref_db
        .iter()
        .filter(|(_, _, code_full, _)| code_full == best_full)
        .count();

    // Si count_ref >= reference_size, devuelve false
    if count_ref >= reference_size {
        return false;
    }

    // 2. Contar en distances cuántas coincidencias apuntan a samples con ese code_full
    let count_dist = distances
        .iter()
        .filter(|(_, target_sample, _)| {
            ref_db
                .iter()
                .find(|(sample, _, code_full, _)| sample == target_sample && code_full == best_full)
                .is_some()
        })
        .count();

    // 3. Calcular proporción y comparar con threshold
    if count_ref == 0 {
        return false;
    }
    let ratio = (count_dist as f64) / (count_ref as f64);
    ratio >= click_threshold
}

/// Busca cliques en niveles taxonómicos específicos y actualiza el objeto Query.
pub fn look_for_cliques(
    conn: &mut Connection,
    distances: &[(String, String, f64)],
    level: usize,
    levels: &[(String, f64)],
    click_size: usize,
    query: &mut Query,
) -> Result<(), duckdb::Error> {
    let clique_search_start = Instant::now();
    
    // 1. Insertar edges
    let edges_insert_start = Instant::now();
    let edges_to_insert: Vec<(String, String, f64)> = distances
        .iter()
        .map(|(q, r, d)| (q.clone(), r.clone(), *d))
        .collect();
    insert_edges(conn, &edges_to_insert)?;
    println!("        ✓ Insertados {} edges [{:.2}ms]", 
             edges_to_insert.len(), edges_insert_start.elapsed().as_millis());

    // 2. IDs únicos
    let mut unique_ids = hashbrown::HashSet::new();
    for (q, r, _) in distances {
        unique_ids.insert(q.clone());
        unique_ids.insert(r.clone());
    }
    let node_ids: Vec<String> = unique_ids.into_iter().collect();

    // 3. Bucle de niveles
    for j in level..levels.len() {
        let level_clique_start = Instant::now();
        let (_level_name, dist_threshold) = &levels[j];
        
        println!(
            "        🔍 Nivel {} umbral {:.4} tamaño mínimo {}",
            j, dist_threshold, click_size
        );

        let clique_check_start = Instant::now();
        if let Some((clique_nodes, clique_edges)) = exists_clique(
            conn,
            &node_ids,
            *dist_threshold,
            click_size,
        )? {
            println!(
                "        ✅ Clique encontrado en nivel {} • nodos: {} • edges: {} [{:.2}ms]",
                j, clique_nodes.len(), clique_edges.len(), clique_check_start.elapsed().as_millis()
            );

            // 4. Actualizar código para el nivel j
            let code_update_start = Instant::now();
            set_clique_code_for_query(conn, query, j)?;
            println!("        ✓ Código de query actualizado [{:.2}ms]", 
                     code_update_start.elapsed().as_millis());
            
            // 5. Propagar valores a nodos del clique
            let clique_update_start = Instant::now();
            update_clique_samples(conn, &clique_nodes, j, query)?;
            println!("        ✓ Samples del clique actualizados [{:.2}ms]", 
                     clique_update_start.elapsed().as_millis());
            
            // 6. Propagar valores a nodos fuera del clique con estado "S"
            let non_clique_update_start = Instant::now();
            update_non_clique_samples(conn, &node_ids, &clique_nodes, j, query)?;
            println!("        ✓ Samples fuera del clique actualizados [{:.2}ms]", 
                     non_clique_update_start.elapsed().as_millis());

            // 7. Eliminar edges inválidos
            let edge_cleanup_start = Instant::now();
            let delete_dist_threshold = if j == levels.len() - 1 {
                1.0
            } else {
                levels[j + 1].1
            };
            
            delete_edges_by_ids(conn, &clique_edges, delete_dist_threshold)?;
            println!(
                "        🗑️ Eliminados {} edges no válidos (dist < {:.4}) [{:.2}ms]",
                clique_edges.len(), delete_dist_threshold, edge_cleanup_start.elapsed().as_millis()
            );

            println!("      ✓ Nivel {} de clique completado [{:.2}ms]", 
                     j, level_clique_start.elapsed().as_millis());
            
            let total_clique_time = clique_search_start.elapsed();
            println!("    ✓ Búsqueda de cliques exitosa [{:.2}ms]", total_clique_time.as_millis());
            return Ok(());
        } else {
            println!("        ○ No hay clique en nivel {} [{:.2}ms]", 
                     j, clique_check_start.elapsed().as_millis());
            break;
        }
    }

    let total_clique_time = clique_search_start.elapsed();
    println!("    ○ Sin cliques desde nivel {} [{:.2}ms]", level, total_clique_time.as_millis());
    Ok(())
}

/// Actualiza el campo code, code_full y code_state de `query` en el nivel `level`
/// basándose en datos de la tabla `code` fusionada.
pub fn set_clique_code_for_query(
    conn: &mut Connection,
    query: &mut Query,
    level: usize,
) -> Result<(), duckdb::Error> {
    let level_int_col = format!("L_{}_int", level);

    // 1. Obtener el valor máximo de code desde la tabla `code`
    let max_code_opt: Option<i64> = if level == 0 {
        conn.query_row(
            &format!("SELECT MAX({}) FROM code", level_int_col),
            [],
            |r| r.get(0),
        )?
    } else {
        let prev_full_col = format!("L_{}_full", level - 1);
        let prev_full = &query.code_full[level - 1];
        conn.query_row(
            &format!(
                "SELECT MAX({level_int_col}) FROM code WHERE {prev_full_col} = ?"
            ),
            params![prev_full],
            |r| r.get(0),
        )?
    };

    // 2. Asignar code = max_code + 1 y estado
    let code_value = max_code_opt.unwrap_or(0) as usize + 1;
    query.code[level] = code_value;
    query.code_state[level] = "C".to_string();

    // 3. Construir code_full
    if level == 0 {
        query.code_full[0] = code_value.to_string();
    } else if code_value == 0 {
        query.code_full[level].clear();
    } else {
        let prev_full = &query.code_full[level - 1];
        if prev_full.is_empty() {
            query.code_full[level].clear();
        } else {
            query.code_full[level] = format!("{}.{}", prev_full, code_value);
        }
    }

    println!(
        "          📝 Query nivel {} actualizado • Code: {} • Full: '{}' • State: C",
        level, query.code[level], query.code_full[level]
    );

    Ok(())
}

/// Propaga los valores de `query` en el nivel `level` a todos los `samples` del clique
/// en la tabla `code` fusionada.
pub fn update_clique_samples(
    conn: &mut Connection,
    clique_nodes: &[String],
    level: usize,
    query: &Query,
) -> Result<(), duckdb::Error> {
    let code_val = query.code[level] as i64;
    let full_val = &query.code_full[level];
    let state_val = &query.code_state[level];
    
    let int_col = format!("L_{}_int", level);
    let full_col = format!("L_{}_full", level);
    let state_col = format!("L_{}_state", level);

    for sample in clique_nodes {
        // Usar UPSERT (INSERT ... ON CONFLICT DO UPDATE)
        conn.execute(
            &format!(
                "INSERT INTO code (sample, {int_col}, {full_col}, {state_col}) \
                 VALUES (?, ?, ?, ?) \
                 ON CONFLICT(sample) DO UPDATE SET \
                 {int_col} = excluded.{int_col}, \
                 {full_col} = excluded.{full_col}, \
                 {state_col} = excluded.{state_col}"
            ),
            params![sample, code_val, full_val, state_val],
        )?;
    }

    Ok(())
}

/// Propaga los valores de `query` en el nivel `level` a los nodos NO pertenecientes al clique
/// con estado "S" en la tabla `code` fusionada.
pub fn update_non_clique_samples(
    conn: &mut Connection,
    all_nodes: &[String],
    clique_nodes: &[String],
    level: usize,
    query: &Query,
) -> Result<(), duckdb::Error> {
    let code_val = query.code[level] as i64;
    let full_val = &query.code_full[level];
    let state_val = "S";
    
    let int_col = format!("L_{}_int", level);
    let full_col = format!("L_{}_full", level);
    let state_col = format!("L_{}_state", level);

    for sample in all_nodes.iter().filter(|s| !clique_nodes.contains(s)) {
        // Usar UPSERT
        conn.execute(
            &format!(
                "INSERT INTO code (sample, {int_col}, {full_col}, {state_col}) \
                 VALUES (?, ?, ?, ?) \
                 ON CONFLICT(sample) DO UPDATE SET \
                 {int_col} = excluded.{int_col}, \
                 {full_col} = excluded.{full_col}, \
                 {state_col} = excluded.{state_col}"
            ),
            params![sample, code_val, full_val, state_val],
        )?;
    }

    Ok(())
}

/// Guarda todos los valores de query en la tabla code fusionada
pub fn update_duckdb(
    conn: &mut Connection,
    query: &Query,
) -> DuckResult<()> {
    let sample = &query.sample_name;
    let n = query.code.len();

    // Construir lista de columnas dinámicamente
    let mut int_cols = Vec::new();
    let mut full_cols = Vec::new();
    let mut state_cols = Vec::new();
    let mut int_values = Vec::new();
    let mut full_values = Vec::new();
    let mut state_values = Vec::new();

    for i in 0..n {
        int_cols.push(format!("L_{}_int", i));
        full_cols.push(format!("L_{}_full", i));
        state_cols.push(format!("L_{}_state", i));
        int_values.push(query.code[i] as i64);
        full_values.push(&query.code_full[i]);
        state_values.push(&query.code_state[i]);
    }

    // Construir SQL dinámicamente
    let all_cols = [&int_cols[..], &full_cols[..], &state_cols[..]].concat();
    let placeholders: String = std::iter::repeat("?").take(all_cols.len() + 1).collect::<Vec<_>>().join(", ");
    
    let update_sets: Vec<String> = all_cols.iter()
        .map(|col| format!("{} = excluded.{}", col, col))
        .collect();

    let sql = format!(
        "INSERT INTO code (sample, {}) VALUES ({}) \
         ON CONFLICT(sample) DO UPDATE SET {}",
        all_cols.join(", "),
        placeholders,
        update_sets.join(", ")
    );

    // Preparar parámetros
    let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();
    params.push(Box::new(sample.clone()));
    
    // Agregar valores en orden: int, full, state
    for val in int_values {
        params.push(Box::new(val));
    }
    for val in full_values {
        params.push(Box::new(val.clone()));
    }
    for val in state_values {
        params.push(Box::new(val.clone()));
    }

    // Ejecutar con params_from_iter
    let mut stmt = conn.prepare(&sql)?;
    stmt.execute(duckdb::params_from_iter(params))?;

    Ok(())
}
