use clap::Args;
use anyhow::{Result, Context, anyhow, bail};
use hashbrown::HashMap;
use std::fs;
use std::path::Path;
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
    
    println!("✓ Cargados {} niveles desde la tabla levels", levels.len());
    Ok(levels)
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
            // Leer como Value genérico y convertir a String
            let val: Option<duckdb::types::Value> = r.get(i)?;
            let val_str = match val {
                Some(duckdb::types::Value::Int(i)) => i.to_string(),
                Some(duckdb::types::Value::Double(f)) => f.to_string(),  // <-- Cambio aquí
                Some(duckdb::types::Value::Text(s)) => s,
                Some(duckdb::types::Value::Boolean(b)) => b.to_string(),
                Some(v) => format!("{:?}", v), // Para otros tipos, usar debug format
                None => String::new(),
            };
            map.insert(col.clone(), val_str);
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

    println!("CPUs: {}", args.cpus);
    println!("Lista de archivos: {}", args.files);

    // Verificar que los archivos existen
    if !Path::new(&args.db).exists() {
        return Err(anyhow::anyhow!("Base de datos no encontrada: {}", args.db));
    }
    if !Path::new(&args.files).exists() {
        return Err(anyhow::anyhow!("Lista de archivos no encontrada: {}", args.files));
    }

    // Abrir conexión mutable a la base de datos
    let mut conn = Connection::open(&args.db)
        .with_context(|| format!("Error abriendo base de datos: {}", args.db))?;

    println!("✓ Conexión a base de datos establecida");

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
    let mut sketch_manager_result = load_sketch_manager_from_db(
        ctx.conn,
        ctx.kmer_size(),
        ctx.sketch_size()
    );

    match &mut sketch_manager_result {
        Ok(sketch_manager) => {
            println!(
                "✓ SketchManager cargado desde DB. Contiene {} sketches",
                sketch_manager.length()
            );

            for line in fs::read_to_string(&args.files)
                .with_context(|| format!("Error leyendo archivo de lista de archivos: {}", args.files))?
                .lines()
            {
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
                    continue;
                }

                update_single_file(fasta_path, &mut ctx, sketch_manager, args.debug)?; // Mut borrow aquí
            }
        }
        Err(ref e) => {
            // Manejo de error; puedes personalizar según tu flujo
            eprintln!("Error al cargar SketchManager desde la base de datos: {e}");
            // return Err(e.clone().into()); // si quieres propagar el error
        }
}



    
    // TODO: Aquí se agregará la lógica para:
    // 1. Comparar sketches query vs reference
    // 2. Actualizar contadores en level_counts
    // 3. Insertar/actualizar registros en tablas code, code_full, code_state
    // 4. Guardar SketchManager actualizado

    println!("\n=== Comando update completado exitosamente ===");
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
    ///
    /// # Arguments
    /// * `path` - Path al archivo FASTA de entrada.
    /// * `ctx` - Contexto que debe poder (por ejemplo) saber el número de niveles.
    pub fn new(path: &Path, ctx: &UpdateCtx) -> anyhow::Result<Self> {
        let sample_name = path.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
        let n_levels = ctx.levels.len();

        // Vector de enteros positivos (inicializados a 0)
        let code = vec![0_usize; n_levels];
        // code_full inicializado con strings vacíos ("")
        let code_full = vec!["".to_string(); n_levels];
        // code_state inicializado con "S" (sin clasificar por defecto)
        let code_state = vec!["S".to_string(); n_levels];
        // Instanciar SketchManager y cargar Sketch
        let mut sketch_manager = SketchManager::new(ctx.kmer_size(), ctx.sketch_size());
        let sketch = Sketch::new(path, ctx.kmer_size(), ctx.sketch_size())?;
        sketch_manager.add_sketch(sketch)?;

        Ok(Query {
            sample_name,
            code,
            code_full,
            code_state,
            sketch: sketch_manager,
        })
    }
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
    _sketch_manager: &mut SketchManager,
    debug: bool
) -> Result<()> {
    // -------- 1. Crear Query (incluye validaciones, sketch y vectores) --------
    let mut query = Query::new(fasta_path, ctx)?;

    println!(
        "✓ Query creado • muestra: {} • k: {} • sketch: {} • niveles: {}",
        query.sample_name,
        ctx.kmer_size(),
        ctx.sketch_size(),
        ctx.num_levels()
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
        // 3.1 Recuperar clasificadores iniciales (condition = "C")
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "C")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "C")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
       // println!(
       //    "Número de referencias para {} en nivel {}: {}",
       //     query.sample_name, i, ref_ids.len()
       // );

        // 3.2 Calcular distancias
        distances = pw_one_to_many(
            &query.sketch,
            _sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        // --- DEBUG: guardar distances si debug == true ---
        if debug {
            let tx = ctx.connection_mut().transaction()?;
            for (q, r, d) in &distances {
                tx.execute(
                    "INSERT INTO debug (Source, Target,dist) VALUES (?1, ?2, ?3)",
                    params![q, r, d],
                )?;
            }
            tx.commit()?;
        }
        // 3.3 Procesar distancias
        if !distances.is_empty() {
            // 5. Actualizar según best_hit
            bh = best_hit(&distances, &ref_db);

            // Obtener valores de best_hit
            let code_val = bh.as_ref().map_or(0, |(_, code, _, _)| *code);
            let code_full_val = bh
                .as_ref()
                .map_or_else(|| "".to_string(), |(_, _, cf, _)| cf.clone());

            if code_val != 0 {
                // Actualizar query con best hit válido
                query.code[i] = code_val;
                query.code_full[i] = code_full_val.clone();

                if is_classifier(
                    &distances,
                    bh.as_ref().unwrap(),
                    &ref_db,
                    click_threshold,
                    reference_size,
                ) {
                    query.code_state[i] = "C".to_string();
                }
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
            "No hay candidatos válidos en nivel {}, buscando con condition = 'ALL'...",
            i
        );
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "ALL")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "ALL")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
        println!("Número de referencias (ALL): {}", ref_ids.len());

        distances = pw_one_to_many(
            &query.sketch,
            _sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        // --- DEBUG: guardar distances si debug == true ---
        if debug {
            let tx = ctx.connection_mut().transaction()?;
            for (q, r, d) in &distances {
                tx.execute(
                    "INSERT INTO debug (Source, Target,dist) VALUES (?1, ?2, ?3)",
                    params![q, r, d],
                )?;
            }
            tx.commit()?;
        }
        println!(
            "✓ Distancias calculadas en ALL con límite {} para {} muestras",
            ctx.levels[i].1,
            distances.len()
        );

        if !distances.is_empty() {
            // 3.5 Buscar cliques y crear nuevos grupos
            let levels = ctx.levels.clone();
            let click_size = ctx.click_size();
            let conn = ctx.connection_mut();
            look_for_cliques(conn, &distances, i, &levels, click_size, &mut query)?;
            break;
        } else {
            break;
        }
    }

    // -------- 4. Actualizar la base de datos y guardar sketch --------
    update_duckdb(ctx.connection_mut(), &query)?;
    let sketch = query
        .sketch
        .get_sketch(&query.sample_name)
        .expect("El sketch debería existir después de Query::new()");
    insert_sketch_object(ctx.conn, &query.sample_name, sketch)
        .with_context(|| format!("Error guardando sketch para muestra: {}", query.sample_name))?;
    _sketch_manager
        .add_sketch(sketch.clone())
        .with_context(|| format!("Error añadiendo sketch para muestra: {}", query.sample_name))?;

    println!(
        "✓ Procesamiento completo • muestra: {} • código final: {:?}",
        query.sample_name, query.code
    );

    Ok(())
}



/// Recupera clasificadores que cumplen condiciones específicas de nivel y grupo.
/// 
/// # Argumentos
/// * `conn` - Conexión a la base de datos DuckDB
/// * `level` - Nivel taxonómico (0, 1, 2, 3, 4...)
/// * `group` - Grupo taxonómico a filtrar (ignorado si level == 0)
/// 
/// # Retorna
/// Vector de tuplas (sample, code, code_full, code_state) donde:
/// - code_state.L_{level} == "C" 
/// - code_full.L_{level} == group (excepto si level == 0)
/// 
/// # Ejemplo
/// ```
/// // Obtener todos los clasificadores de nivel 0 (ignora group)
/// let results = retrieve_classifiers(conn, 0, "0")?;
/// 
/// // Obtener clasificadores de nivel 2 del grupo "Escherichia"
/// let results = retrieve_classifiers(conn, 2, "1.1")?;
/// ```

pub fn retrieve_classifiers(
    conn: &Connection,
    level: usize,
    group: &str,
    condition: &str,
) -> Result<Vec<(String, usize, String, String)>> {
    
    let level_col = format!("L_{}", level);
    
    // Para los filtros, necesitamos usar level-1 (excepto cuando level = 0)
    let filter_level_col = if level == 0 {
        format!("L_{}", level)  // Si level es 0, usamos L_0 para ambos filtros
    } else {
        format!("L_{}", level - 1)  // Si level > 0, usamos L_(level-1) para ambos filtros
    };
    
    let (sql, params): (String, Vec<&dyn duckdb::ToSql>) = match (level, condition) {
        // Level 0 con condición 'C'
        (0, "C") => (
            format!(
                "SELECT 
                    c.sample,
                    c.{level_col} as code,
                    cf.{level_col} as code_full,
                    cs.{level_col} as code_state
                FROM code c
                JOIN code_full cf ON c.sample = cf.sample  
                JOIN code_state cs ON c.sample = cs.sample
                WHERE cs.{level_col} = 'C'"  // ← Cambio: cs.{level_col} en lugar de cs.{filter_level_col}
            ),
            Vec::new()
        ),
        
        // Level 0 con condición 'ALL'
        (0, "ALL") => (
            format!(
                "SELECT 
                    c.sample,
                    c.{level_col} as code,
                    cf.{level_col} as code_full,
                    cs.{level_col} as code_state
                FROM code c
                JOIN code_full cf ON c.sample = cf.sample  
                JOIN code_state cs ON c.sample = cs.sample"
            ),
            Vec::new()
        ),
        
        // Level > 0 con condición 'C'
        (level, "C") if level > 0 => (
            format!(
                "SELECT 
                    c.sample,
                    c.{level_col} as code,
                    cf.{level_col} as code_full,
                    cs.{level_col} as code_state
                FROM code c
                JOIN code_full cf ON c.sample = cf.sample  
                JOIN code_state cs ON c.sample = cs.sample
                WHERE cs.{level_col} = 'C' 
                AND cf.{filter_level_col} = ?"  // ← cs.{level_col} + cf.{filter_level_col}
            ),
            vec![&group as &dyn duckdb::ToSql]
        ),
        
        // Level > 0 con condición 'ALL'
        (level, "ALL") if level > 0 => (
            format!(
                "SELECT 
                    c.sample,
                    c.{level_col} as code,
                    cf.{level_col} as code_full,
                    cs.{level_col} as code_state
                FROM code c
                JOIN code_full cf ON c.sample = cf.sample  
                JOIN code_state cs ON c.sample = cs.sample
                WHERE cf.{filter_level_col} = ?"  // ← Solo cf.{filter_level_col}
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
        //println!("Processing row {}", row_count);
        
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
                "UNKNOWN".to_string()  // ← Valor por defecto para NULL
            }
        };
        
        let code_state: String = match row.get::<_, Option<String>>(3)? {
            Some(val) => val,
            None => {
                println!("WARNING: NULL code_state at row {}, using 'UNKNOWN'", row_count);
                "UNKNOWN".to_string()  // ← Valor por defecto para NULL
            }
        };
        
        
        result.push((sample, code, code_full, code_state));
    }

    Ok(result)
}



/// Encuentra la similitud maxima y devuelve la tupla correspondiente de ref_db.
/// 
/// # Argumentos
/// * `distances` - Vector de tuplas (query_name, ref_name, distance)
/// * `ref_db` - Vector de tuplas (sample, code, code_full, code_state) de retrieve_classifiers
/// 
/// # Retorna
/// Option con la tupla de ref_db que tiene el sample correspondiente al mejor hit,
/// o None si no se encuentra coincidencia o no hay distancias válidas.
pub fn best_hit(
    distances: &[(String, String, f64)],
    ref_db: &[(String, usize, String, String)],
) -> Option<(String, usize, String, String)> {
    // 1. Encontrar la distancia mínima (ignorar NaN)
    let best_distance = distances
        .iter()
        .filter(|(_, _, dist)| !dist.is_nan())  // Filtrar NaN
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())?;  // Comparar por distancia (índice 2)
    //println!("Mejor hit encontrado: {} con distancia {:.4}", best_distance.1, best_distance.2);
    // 2. Obtener el ref_name del mejor hit
    let best_ref_name = &best_distance.1;  // Segunda columna de distances
    
    // 3. Buscar la tupla correspondiente en ref_db
    ref_db
        .iter()
        .find(|(sample, _, _, _)| sample == best_ref_name)  // Primera columna de ref_db
        .cloned()  // Clonar para devolver la tupla completa
}

/// Determina si un best_hit pasa el criterio de clasificación.
///
/// # Argumentos
/// * `distances`       - Vector de tuplas (query_name, ref_name, dist).
/// * `best_hit`        - Tupla (sample, code, code_full, code_state) del mejor hit.
/// * `ref_db`          - Vector de tuplas (sample, code, code_full, code_state).
/// * `click_threshold` - Umbral mínimo de proporción (p.ej. 0.025).
/// * `reference_size`  - Tamaño máximo permitido para el grupo; si se iguala o supera, falla.
///
/// # Retorna
/// `true` si (count_dist / count_ref) >= click_threshold **y** count_ref < reference_size;
/// `false` en caso contrario.
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
        // Evitar división por cero; si no hay references, falla
        return false;
    }
    let ratio = (count_dist as f64) / (count_ref as f64);
    ratio >= click_threshold
}
/// Busca cliques en niveles taxonómicos específicos y actualiza el objeto Query.
///
/// # Argumentos
/// * `conn` - Conexión mutable a DuckDB
/// * `distances` - Vector de distancias (query_name, ref_name, distance) a insertar
/// * `level` - Nivel inicial desde donde empezar la búsqueda (inclusive)
/// * `ctx` - Contexto con información de niveles y parámetros
/// * `query` - Objeto Query que será modificado con los resultados del clique
///
/// # Proceso
/// 1. Inserta las distancias en la tabla edges
/// 2. Busca cliques desde el nivel especificado hacia arriba
/// 3. Si encuentra un clique, actualiza el Query y elimina edges
///
/// # Retorna
/// Result con () en caso de éxito, o Error en caso de fallo
pub fn look_for_cliques(
    conn: &mut Connection,
    distances: &[(String, String, f64)],
    level: usize,
    levels: &[(String, f64)],
    click_size: usize,
    query: &mut Query,
) -> Result<(), duckdb::Error> {
    // 1. Insertar edges
    let edges_to_insert: Vec<(String, String, f64)> = distances
        .iter()
        .map(|(q, r, d)| (q.clone(), r.clone(), *d))
        .collect();
    insert_edges(conn, &edges_to_insert)?;
    println!("✓ Insertados {} edges", edges_to_insert.len());

    // 2. IDs únicos
    let mut unique_ids = hashbrown::HashSet::new();
    for (q, r, _) in distances {
        unique_ids.insert(q.clone());
        unique_ids.insert(r.clone());
    }
    let node_ids: Vec<String> = unique_ids.into_iter().collect();

    // 3. Bucle de niveles
    for j in level..levels.len() {
        let (_level_name, dist_threshold) = &levels[j];
        
        println!(
            "🔍 Nivel {} umbral {:.4} tamaño mínimo {}",
            j, dist_threshold, click_size
        );

        if let Some((clique_nodes, clique_edges)) = exists_clique(
            conn,
            &node_ids,
            *dist_threshold,
            click_size,
        )? {
            println!(
                "✅ Clique en nivel {} • nodos: {} • edges: {}",
                j,
                clique_nodes.len(),
                clique_edges.len()
            );

            // 4. Actualizar código para el nivel j
            set_clique_code_for_query(conn, query, j)?;
            // 5. Propagar valores a nodos del clique
            update_clique_samples(conn, &clique_nodes, j, query)?;
            // 6. Propagar valores a nodos fuera del clique con estado "S"
            update_non_clique_samples(conn, &node_ids, &clique_nodes, j, query)?;

            // 7. Eliminar edges inválidos: solo aquellos del clique con dist < dist_threshold
            let delete_dist_threshold = if j == levels.len() - 1 {
                // Si es el último nivel, usar el threshold del siguiente nivel (que no existe)
                1.0 // Valor arbitrario alto para eliminar todos los edges
            } else {
                levels[j + 1].1 // Umbral del siguiente nivel
            };
            
            delete_edges_by_ids(conn, &clique_edges, delete_dist_threshold)?;
            println!(
                "🗑️ Eliminados {} edges no válidos (dist < {:.4})",
                clique_edges.len(),
                delete_dist_threshold
            );

            return Ok(());
        } else {
            println!("○ No hay clique en nivel {}", j);
            break; // No hay clique, salir del bucle
        }
    }

    println!("○ Sin cliques desde nivel {}", level);
    Ok(())
}




/// Actualiza el campo code, code_full y code_state de `query` en el nivel `level`
/// basándose en datos de la tabla `code`.
///
/// - Si `level == 0`, obtiene el máximo `code.L_0` de la tabla `code`.
/// - Si `level > 0`, obtiene el máximo `code.L_{level}`
///   para filas cuyo `code_full.L_{level-1}` coincida con `query.code_full[level-1]`.
/// - Asigna `query.code[level]` y `query.code_full[level]` directamente.
/// - Marca `query.code_state[level] = "C"`.
pub fn set_clique_code_for_query(
    conn: &mut Connection,
    query: &mut Query,
    level: usize,
) -> Result<(), duckdb::Error> {
    let level_col = format!("L_{}", level);

    // 1. Obtener el valor máximo de code desde la tabla `code`
    let max_code_opt: Option<i64> = if level == 0 {
        conn.query_row(
            &format!("SELECT MAX({}) FROM code", level_col),
            [],
            |r| r.get(0),
        )?
    } else {
        let prev_col = format!("L_{}", level - 1);
        let prev_full = &query.code_full[level - 1];
        conn.query_row(
            &format!(
                "SELECT MAX(c.{level_col}) \
                 FROM code c \
                 JOIN code_full cf ON c.sample = cf.sample \
                 WHERE cf.{prev_col} = ?"
            ),
            params![prev_full],
            |r| r.get(0),
        )?
    };

    // 2. Asignar code = max_code + 1 (o 1 si no hay resultado) y estado
    let code_value = max_code_opt.unwrap_or(0) as usize + 1;
    query.code[level] = code_value;
    query.code_state[level] = "C".to_string();

    // 3. Construir code_full inline
    if level == 0 {
        query.code_full[0] = code_value.to_string();
    } else if code_value == 0 {
        // Nunca sucede porque code_value >= 1
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
        "📝 Query nivel {} actualizado • Code: {} • Full: '{}' • State: C",
        level, query.code[level], query.code_full[level]
    );

    Ok(())
}

/// Propaga los valores de `query` en el nivel `level` a todos los `samples` del clique,
/// actualizando e insertando en las tablas `code`, `code_full` y `code_state`.
///
/// # Argumentos
/// * `conn`         – Conexión mutable a DuckDB.
/// * `clique_nodes` – Slice de IDs de muestras a actualizar.
/// * `level`        – Nivel de taxonomía (ej. 0,1,2…).
/// * `query`        – Objeto Query con los valores ya asignados para ese nivel.
///
/// # Retorna
/// Result<(), duckdb::Error>
pub fn update_clique_samples(
    conn: &mut Connection,
    clique_nodes: &[String],
    level: usize,
    query: &Query,
) -> Result<(), duckdb::Error> {
    let code_val = query.code[level] as i64;
    let full_val = &query.code_full[level];
    let state_val = &query.code_state[level];
    let level_col   = format!("L_{}", level);

    for sample in clique_nodes {
        // UPDATE existentes
        conn.execute(
            &format!("UPDATE code SET {} = ? WHERE sample = ?", level_col),
            params![code_val, sample],
        )?;
        conn.execute(
            &format!("UPDATE code_full SET {} = ? WHERE sample = ?", level_col),
            params![full_val, sample],
        )?;
        conn.execute(
            &format!("UPDATE code_state SET {} = ? WHERE sample = ?", level_col),
            params![state_val, sample],
        )?;

        // INSERT si no existen
        conn.execute(
            &format!(
                "INSERT INTO code (sample, {}) \
                 SELECT ?, ? WHERE NOT EXISTS (SELECT 1 FROM code WHERE sample = ?)",
                level_col
            ),
            params![sample, code_val, sample],
        )?;
        conn.execute(
            &format!(
                "INSERT INTO code_full (sample, {}) \
                 SELECT ?, ? WHERE NOT EXISTS (SELECT 1 FROM code_full WHERE sample = ?)",
                level_col
            ),
            params![sample, full_val, sample],
        )?;
        conn.execute(
            &format!(
                "INSERT INTO code_state (sample, {}) \
                 SELECT ?, ? WHERE NOT EXISTS (SELECT 1 FROM code_state WHERE sample = ?)",
                level_col
            ),
            params![sample, state_val, sample],
        )?;
    }

    Ok(())
}

/// Propaga los valores de `query` en el nivel `level` a los nodos NO pertenecientes al clique,
/// actualizando e insertando en las tablas `code`, `code_full` y `code_state` con estado "S".
///
/// # Argumentos
/// * `conn`       – Conexión mutable a DuckDB.
/// * `all_nodes`  – Slice de todos los IDs de muestras consideradas.
/// * `clique_nodes` – Slice de IDs de muestras que forman el clique.
/// * `level`      – Nivel de taxonomía (ej. 0,1,2…).
/// * `query`      – Objeto Query con los valores ya asignados para ese nivel.
///
/// # Retorna
/// Result<(), duckdb::Error>
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
    let level_col   = format!("L_{}", level);

    for sample in all_nodes.iter().filter(|s| !clique_nodes.contains(s)) {
        // UPDATE existentes
        conn.execute(
            &format!("UPDATE code SET {} = ? WHERE sample = ?", level_col),
            params![code_val, sample],
        )?;
        conn.execute(
            &format!("UPDATE code_full SET {} = ? WHERE sample = ?", level_col),
            params![full_val, sample],
        )?;
        conn.execute(
            &format!("UPDATE code_state SET {} = ? WHERE sample = ?", level_col),
            params![state_val, sample],
        )?;

        // INSERT si no existen
        conn.execute(
            &format!(
                "INSERT INTO code (sample, {}) \
                 SELECT ?, ? WHERE NOT EXISTS (SELECT 1 FROM code WHERE sample = ?)",
                level_col
            ),
            params![sample, code_val, sample],
        )?;
        conn.execute(
            &format!(
                "INSERT INTO code_full (sample, {}) \
                 SELECT ?, ? WHERE NOT EXISTS (SELECT 1 FROM code_full WHERE sample = ?)",
                level_col
            ),
            params![sample, full_val, sample],
        )?;
        conn.execute(
            &format!(
                "INSERT INTO code_state (sample, {}) \
                 SELECT ?, ? WHERE NOT EXISTS (SELECT 1 FROM code_state WHERE sample = ?)",
                level_col
            ),
            params![sample, state_val, sample],
        )?;
    }

    Ok(())
}



/// Guarda todos los valores de query (code, code_full, code_state)
/// en las tablas code, code_full y code_state de DuckDB para el sample actual.
///
/// Requiere que Query tenga los campos:
/// - sample_name: String
/// - code: Vec<usize>
/// - code_full: Vec<String>
/// - code_state: Vec<String>
pub fn update_duckdb(
    conn: &mut Connection,
    query: &Query,
) -> DuckResult<()> {
    // Asume que los vectores code, code_full y code_state tienen igual longitud
    let sample = &query.sample_name;
    let n = query.code.len();

    for i in 0..n {
        let col = format!("L_{}", i);

        // GUARDAR EN TABLA "code"
        conn.execute(
            &format!(
                "INSERT INTO code (sample, {col}) VALUES (?, ?) \
                 ON CONFLICT(sample) DO UPDATE SET {col} = excluded.{col}"
            ),
            params![sample, query.code[i] as i64], // DuckDB espera i64 para enteros
        )?;

        // GUARDAR EN TABLA "code_full"
        conn.execute(
            &format!(
                "INSERT INTO code_full (sample, {col}) VALUES (?, ?) \
                 ON CONFLICT(sample) DO UPDATE SET {col} = excluded.{col}"
            ),
            params![sample, &query.code_full[i]],
        )?;

        // GUARDAR EN TABLA "code_state"
        conn.execute(
            &format!(
                "INSERT INTO code_state (sample, {col}) VALUES (?, ?) \
                 ON CONFLICT(sample) DO UPDATE SET {col} = excluded.{col}"
            ),
            params![sample, &query.code_state[i]],
        )?;
    }

    Ok(())
}
