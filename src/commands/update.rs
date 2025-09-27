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

/// Arguments for the update subcommand
#[derive(Args, Debug)]
pub struct UpdateArgs {
    /// Path to DuckDB database
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,

    /// Number of CPUs to parallelize with rayon
    #[arg(long, value_name = "N_CPUS", default_value_t = 1)]
    pub cpus: usize,

    /// Plain text file with paths to FASTA files (one per line)
    #[arg(long, required = true, value_name = "FILES_LIST")]
    pub files: String,

    #[arg(long, required= false, value_name = "DEBUG")]
    pub debug: bool
}

/// Context with preloaded information needed for the whole `update` workflow
pub struct UpdateCtx<'a> {
    /// Mutable connection to DuckDB (reusable in all phases)
    pub conn: &'a mut Connection,
    /// Levels table (L_0 ... L_N → dist) as ordered vector
    pub levels: Vec<(String, f64)>,
    /// Entire metadata row as (column → value in UTF-8)
    pub metadata_map: HashMap<String, String>,
}

impl<'a> UpdateCtx<'a> {
    /// Convenience accessors for metadata
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

    /// Access to mutable connection
    pub fn connection_mut(&mut self) -> &mut Connection {
        self.conn
    }

    pub fn reference_size(&self) -> usize {
        self.metadata_map["reference_size"].parse().unwrap_or(100)
    }

    /// Number of levels
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    /// Get distance for a specific level
    pub fn get_level_distance(&self, level: &str) -> Option<f64> {
        self.levels.iter().find(|(l, _)| l == level).map(|(_, d)| *d)
    }
}

/// Reads the `levels` table and returns a Vec<(String, f64)> ordered by level
fn load_levels_vec(conn: &Connection) -> Result<Vec<(String, f64)>> {
    let start = Instant::now();

    let mut stmt = conn.prepare("SELECT level, dist FROM levels ORDER BY dist ASC")
        .context("Error preparing query for levels table")?;

    let mut rows = stmt.query([])
        .context("Error running query for levels table")?;

    let mut levels = Vec::new();
    while let Some(row) = rows.next()? {
        let level: String = row.get(0)
            .context("Error getting level column")?;
        let dist: f64 = row.get(1)
            .context("Error getting dist column")?;
        levels.push((level, dist));
    }

    let duration = start.elapsed();
    println!("✓ Loaded {} levels from levels table [{:.2}ms]",
             levels.len(), duration.as_millis());
    Ok(levels)
}

/// Reads metadata (single row) → HashMap<column, value_as_text>
fn load_metadata_map(conn: &Connection) -> Result<HashMap<String, String>> {
    let start = Instant::now();

    // Get columns dynamically
    let mut cols_stmt = conn.prepare(
        "SELECT column_name FROM information_schema.columns \
         WHERE table_name = 'metadata' ORDER BY ordinal_position"
    ).context("Error preparing metadata columns query")?;

    let mut cols_rows = cols_stmt.query([])
        .context("Error running metadata columns query")?;

    let mut cols = Vec::new();
    while let Some(row) = cols_rows.next()? {
        let col_name: String = row.get(0)
            .context("Error getting column name")?;
        cols.push(col_name);
    }

    // Build dynamic SELECT query
    let select_sql = format!("SELECT {} FROM metadata LIMIT 1", cols.join(", "));
    let mut row_stmt = conn.prepare(&select_sql)
        .context("Error preparing metadata query")?;

    let row = row_stmt.query_row([], |r| {
        let mut map = HashMap::new();
        for (i, col) in cols.iter().enumerate() {
            // Read as generic Value and convert to String
            let val: Option<duckdb::types::Value> = r.get(i)?;
            let val_str = match val {
                Some(duckdb::types::Value::Int(i)) => i.to_string(),
                Some(duckdb::types::Value::Double(f)) => f.to_string(),
                Some(duckdb::types::Value::Text(s)) => s,
                Some(duckdb::types::Value::Boolean(b)) => b.to_string(),
                Some(v) => format!("{:?}", v),
                None => String::new(),
            };
            map.insert(col.clone(), val_str);
        }
        Ok(map)
    }).context("Error reading metadata row")?;

    let duration = start.elapsed();
    println!("✓ Loaded {} fields from metadata table [{:.2}ms]",
             row.len(), duration.as_millis());
    Ok(row)
}

/// Main function for the update command
pub fn update_command(args: &UpdateArgs) -> Result<()> {
    let command_start = Instant::now();

    println!("=== Starting update command ===");
    println!("Database: {}", args.db);
    println!("CPUs: {}", args.cpus);
    println!("File list: {}", args.files);

    // Check that files exist
    let validation_start = Instant::now();
    if !Path::new(&args.db).exists() {
        return Err(anyhow::anyhow!("Database not found: {}", args.db));
    }
    if !Path::new(&args.files).exists() {
        return Err(anyhow::anyhow!("Files list not found: {}", args.files));
    }
    println!("✓ Files validated [{:.2}ms]", validation_start.elapsed().as_millis());

    // Open mutable connection to database
    let db_connect_start = Instant::now();
    let mut conn = Connection::open(&args.db)
        .with_context(|| format!("Error opening database: {}", args.db))?;
    println!("✓ Database connection established [{:.2}ms]",
             db_connect_start.elapsed().as_millis());

    // Preload levels and metadata tables
    let levels = load_levels_vec(&conn)?;
    let metadata_map = load_metadata_map(&conn)?;

    // Create shared context
    let mut ctx = UpdateCtx {
        conn: &mut conn,
        levels,
        metadata_map,
    };

    println!("\n=== Project information ===");
    println!("Genus: {}", ctx.genus());
    println!("Acronym: {}", ctx.acronym());
    println!("K-mer size: {}", ctx.kmer_size());
    println!("Sketch size: {}", ctx.sketch_size());
    println!("Click size: {}", ctx.click_size());
    println!("Click threshold: {}", ctx.click_threshold());
    println!("Available levels: {}", ctx.levels.len());

    // Load SketchManager from DuckDB using UpdateCtx
    let sketch_load_start = Instant::now();
    let mut sketch_manager_result = load_sketch_manager_from_db(
        ctx.conn,
        ctx.kmer_size(),
        ctx.sketch_size()
    );
    println!("✓ SketchManager loaded from DB [{:.2}ms]",
             sketch_load_start.elapsed().as_millis());

    match &mut sketch_manager_result {
        Ok(sketch_manager) => {
            println!(
                "✓ SketchManager loaded from DB. Contains {} sketches",
                sketch_manager.length()
            );

            // Read file list
            let file_read_start = Instant::now();
            let file_content = fs::read_to_string(&args.files)
                .with_context(|| format!("Error reading file list: {}", args.files))?;
            let files: Vec<&str> = file_content.lines().collect();
            println!("✓ File list read: {} files [{:.2}ms]",
                     files.len(), file_read_start.elapsed().as_millis());

            // Process each file
            let mut processed_files = 0;
            let mut skipped_files = 0;
            let processing_start = Instant::now();

            for line in files {
                let fasta_path = Path::new(line.trim());
                if !fasta_path.exists() {
                    eprintln!("FASTA file not found: {}", fasta_path.display());
                    continue;
                }

                let sample_name = fasta_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string();

                if sketch_manager.contains(&sample_name) {
                    println!("✓ Sketch for {} already exists, skipping.", sample_name);
                    skipped_files += 1;
                    continue;
                }

                let file_start = Instant::now();
                update_single_file(fasta_path, &mut ctx, sketch_manager, args.debug)?;
                processed_files += 1;
                println!("  └─ File processed [{:.2}s]", file_start.elapsed().as_secs_f32());
            }

            println!("✓ File processing completed: {} processed, {} skipped [{:.2}s]",
                     processed_files, skipped_files, processing_start.elapsed().as_secs_f32());
        }
        Err(ref e) => {
            eprintln!("Error loading SketchManager from database: {e}");
        }
    }

    let total_duration = command_start.elapsed();
    println!("\n=== Update command completed successfully [{:.2}s] ===",
             total_duration.as_secs_f32());
    Ok(())
}

pub struct Query {
    pub sample_name: String,
    pub code: Vec<usize>,         // vector of positive integers
    pub code_full: Vec<String>,   // vector of strings
    pub code_state: Vec<String>,  // vector of strings
    pub sketch: SketchManager,
}

impl Query {
    /// Creates a new Query, loading SketchManager+Sketch and vectors of appropriate size
    pub fn new(path: &Path, ctx: &UpdateCtx) -> anyhow::Result<Self> {
        let query_creation_start = Instant::now();

        let sample_name = path.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
        let n_levels = ctx.levels.len();

        // Vector of positive integers (initialized to 0)
        let code = vec![0_usize; n_levels];
        // code_full initialized with empty strings ("")
        let code_full = vec!["".to_string(); n_levels];
        // code_state initialized with "S" (unclassified by default)
        let code_state = vec!["S".to_string(); n_levels];

        // Instantiate SketchManager and load Sketch
        let sketch_creation_start = Instant::now();
        let mut sketch_manager = SketchManager::new(ctx.kmer_size(), ctx.sketch_size());
        let sketch = Sketch::new(path, ctx.kmer_size(), ctx.sketch_size())?;
        sketch_manager.add_sketch(sketch)?;
        println!("    ✓ Sketch created [{:.2}ms]", sketch_creation_start.elapsed().as_millis());

        println!("    ✓ Query initialized [{:.2}ms]", query_creation_start.elapsed().as_millis());

        Ok(Query {
            sample_name,
            code,
            code_full,
            code_state,
            sketch: sketch_manager,
        })
    }
}

/// Processes **one** FASTA file
pub fn update_single_file(
    fasta_path: &Path,
    ctx: &mut UpdateCtx,
    _sketch_manager: &mut SketchManager,
    debug: bool
) -> Result<()> {
    let file_processing_start = Instant::now();

    // -------- 1. Create Query (includes validations, sketch and vectors) --------
    let query_start = Instant::now();
    let mut query = Query::new(fasta_path, ctx)?;
    println!("  ✓ Query created • sample: {} • k: {} • sketch: {} • levels: {} [{:.2}ms]",
        query.sample_name,
        ctx.kmer_size(),
        ctx.sketch_size(),
        ctx.num_levels(),
        query_start.elapsed().as_millis()
    );

    // -------- 2. Processing variables --------
    let mut ref_db: Vec<(String, usize, String, String)>;
    let mut distances: Vec<(String, String, f64)>;
    let mut bh: Option<(String, usize, String, String)> = None;
    let mut ref_ids: Vec<String> = Vec::new();

    let click_threshold = ctx.click_threshold();
    let click_size = ctx.click_size();
    let reference_size = ctx.reference_size();

    // -------- 3. Process each level in ascending order of distance --------
    for i in 0..ctx.levels.len() {
        let level_start = Instant::now();
        println!("    🔍 Processing level {} [{:.4} threshold]", i, ctx.levels[i].1);

        // 3.1 Retrieve initial classifiers (condition = "C")
        let classifiers_start = Instant::now();
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "C")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "C")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
        println!("      ✓ Classifiers obtained: {} [{:.2}ms]",
                 ref_ids.len(), classifiers_start.elapsed().as_millis());

        // 3.2 Calculate distances
        let distances_start = Instant::now();
        distances = pw_one_to_many(
            &query.sketch,
            _sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        println!("      ✓ Distances calculated: {} [{:.2}ms]",
                 distances.len(), distances_start.elapsed().as_millis());

        // --- DEBUG: save distances if debug == true ---
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
            println!("      ✓ Debug data saved [{:.2}ms]", debug_start.elapsed().as_millis());
        }

        // 3.3 Process distances
        if !distances.is_empty() {
            let best_hit_start = Instant::now();
            // 5. Update according to best_hit
            bh = best_hit(&distances, &ref_db);
            println!("      ✓ Best hit calculated [{:.2}ms]", best_hit_start.elapsed().as_millis());

            // Get values for best_hit
            let code_val = bh.as_ref().map_or(0, |(_, code, _, _)| *code);
            let code_full_val = bh
                .as_ref()
                .map_or_else(|| "".to_string(), |(_, _, cf, _)| cf.clone());

            if code_val != 0 {
                // Update query with valid best hit
                query.code[i] = code_val;
                query.code_full[i] = code_full_val.clone();

                if i > 0 && query.code_full[i].starts_with(query.code_full[i-1].as_str()) == false {
                    println!("Warning: code_full of best hit '{}' is inconsistent with previous level '{}' at level {}, treating as no candidates",
                        query.code_full[i], query.code_full[i-1], i);
                    panic!("Inconsistency in code_full between consecutive levels");
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
                println!("      ✓ Classifier evaluated [{:.2}ms]", classifier_check_start.elapsed().as_millis());

                println!("    ✓ Level {} completed with best hit [{:.2}ms]", i, level_start.elapsed().as_millis());
                // Continue to next level
                continue;
            } else {
                println!(
                    "Warning: Best hit code = 0 at level {}, treating as no candidates",
                    i
                );
            }
        }

        // 3.4 No valid candidates (or code == 0): search with condition = "ALL"
        println!(
            "      No valid candidates at level {}, searching with condition = 'ALL'...",
            i
        );

        let all_classifiers_start = Instant::now();
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "ALL")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "ALL")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
        println!("      ✓ References (ALL) obtained: {} [{:.2}ms]",
                 ref_ids.len(), all_classifiers_start.elapsed().as_millis());

        let all_distances_start = Instant::now();
        distances = pw_one_to_many(
            &query.sketch,
            _sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        println!("      ✓ Distances (ALL) calculated: {} [{:.2}ms]",
                 distances.len(), all_distances_start.elapsed().as_millis());

        // --- DEBUG: save distances if debug == true ---
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
            println!("      ✓ Debug data (ALL) saved [{:.2}ms]", debug_start.elapsed().as_millis());
        }

        if !distances.is_empty() {
            // 3.5 Search for cliques and create new groups
            let cliques_start = Instant::now();
            let levels = ctx.levels.clone();
            let click_size = ctx.click_size();
            let conn = ctx.connection_mut();
            look_for_cliques(conn, &distances, i, &levels, click_size, &mut query)?;
            println!("      ✓ Clique search completed [{:.2}ms]", cliques_start.elapsed().as_millis());
            println!("    ✓ Level {} completed with cliques [{:.2}ms]", i, level_start.elapsed().as_millis());
            break;
        } else {
            println!("    ○ Level {} without valid distances [{:.2}ms]", i, level_start.elapsed().as_millis());
            break;
        }
    }

    // -------- 4. Update database and save sketch --------
    let db_update_start = Instant::now();
    update_duckdb(ctx.connection_mut(), &query)?;
    println!("    ✓ Database updated [{:.2}ms]", db_update_start.elapsed().as_millis());

    let sketch_save_start = Instant::now();
    let sketch = query
        .sketch
        .get_sketch(&query.sample_name)
        .expect("Sketch should exist after Query::new()");
    insert_sketch_object(ctx.conn, &query.sample_name, sketch)
        .with_context(|| format!("Error saving sketch for sample: {}", query.sample_name))?;
    _sketch_manager
        .add_sketch(sketch.clone())
        .with_context(|| format!("Error adding sketch for sample: {}", query.sample_name))?;
    println!("    ✓ Sketch saved [{:.2}ms]", sketch_save_start.elapsed().as_millis());

    println!(
        "✓ Complete processing • sample: {} • final code: {:?} [{:.2}s]",
        query.sample_name, query.code, file_processing_start.elapsed().as_secs_f32()
    );

    Ok(())
}

/// Retrieves classifiers from the merged code table
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

    // For filters, need level-1 (except when level = 0)
    let filter_full_col = if level == 0 {
        format!("L_{}_full", level)
    } else {
        format!("L_{}_full", level - 1)
    };

    let (sql, params): (String, Vec<&dyn duckdb::ToSql>) = match (level, condition) {
        // Level 0 with condition 'C'
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

        // Level 0 with condition 'ALL'
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

        // Level > 0 with condition 'C'
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

        // Level > 0 with condition 'ALL'
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

        // Unrecognized condition
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

        // Safe handling of nullable columns
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
    if duration.as_millis() > 10 { // Only report if takes more than 10ms
        println!("        └─ SQL query executed [{:.2}ms]", duration.as_millis());
    }

    Ok(result)
}

/// Finds the maximum similarity and returns the corresponding tuple from ref_db.
pub fn best_hit(
    distances: &[(String, String, f64)],
    ref_db: &[(String, usize, String, String)],
) -> Option<(String, usize, String, String)> {
    // 1. Find minimum distance (ignore NaN)
    let best_distance = distances
        .iter()
        .filter(|(_, _, dist)| !dist.is_nan())
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())?;

    // 2. Get ref_name of best hit
    let best_ref_name = &best_distance.1;

    // 3. Find the corresponding tuple in ref_db
    ref_db
        .iter()
        .find(|(sample, _, _, _)| sample == best_ref_name)
        .cloned()
}

/// Determines if a best_hit passes classifier criteria.
pub fn is_classifier(
    distances: &[(String, String, f64)],
    best_hit: &(String, usize, String, String),
    ref_db: &[(String, usize, String, String)],
    click_threshold: f64,
    reference_size: usize,
) -> bool {
    // Extract code_full from best_hit
    let best_full = &best_hit.2;

    // 1. Count how many rows in ref_db have same code_full
    let count_ref = ref_db
        .iter()
        .filter(|(_, _, code_full, _)| code_full == best_full)
        .count();

    // If count_ref >= reference_size, return false
    if count_ref >= reference_size {
        return false;
    }

    // 2. Count in distances how many matches point to samples with that code_full
    let count_dist = distances
        .iter()
        .filter(|(_, target_sample, _)| {
            ref_db
                .iter()
                .find(|(sample, _, code_full, _)| sample == target_sample && code_full == best_full)
                .is_some()
        })
        .count();

    // 3. Compute ratio and compare to threshold
    if count_ref == 0 {
        return false;
    }
    let ratio = (count_dist as f64) / (count_ref as f64);
    ratio >= click_threshold
}

/// Searches for cliques at specific taxonomic levels and updates the Query object.
pub fn look_for_cliques(
    conn: &mut Connection,
    distances: &[(String, String, f64)],
    level: usize,
    levels: &[(String, f64)],
    click_size: usize,
    query: &mut Query,
) -> Result<(), duckdb::Error> {
    let clique_search_start = Instant::now();

    // 1. Insert edges
    let edges_insert_start = Instant::now();
    let edges_to_insert: Vec<(String, String, f64)> = distances
        .iter()
        .map(|(q, r, d)| (q.clone(), r.clone(), *d))
        .collect();
    insert_edges(conn, &edges_to_insert)?;
    println!("        ✓ Inserted {} edges [{:.2}ms]",
             edges_to_insert.len(), edges_insert_start.elapsed().as_millis());

    // 2. Unique IDs
    let mut unique_ids = hashbrown::HashSet::new();
    for (q, r, _) in distances {
        unique_ids.insert(q.clone());
        unique_ids.insert(r.clone());
    }
    let node_ids: Vec<String> = unique_ids.into_iter().collect();

    // 3. Level loop
    for j in level..levels.len() {
        let level_clique_start = Instant::now();
        let (_level_name, dist_threshold) = &levels[j];

        println!(
            "        🔍 Level {} threshold {:.4} minimum size {}",
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
                "        ✅ Clique found at level {} • nodes: {} • edges: {} [{:.2}ms]",
                j, clique_nodes.len(), clique_edges.len(), clique_check_start.elapsed().as_millis()
            );

            // 4. Update code for level j
            let code_update_start = Instant::now();
            set_clique_code_for_query(conn, query, j)?;
            println!("        ✓ Query code updated [{:.2}ms]",
                     code_update_start.elapsed().as_millis());

            // 5. Propagate values to clique nodes
            let clique_update_start = Instant::now();
            update_clique_samples(conn, &clique_nodes, j, query)?;
            println!("        ✓ Clique samples updated [{:.2}ms]",
                     clique_update_start.elapsed().as_millis());

            // 6. Propagate values to nodes outside clique with state "S"
            let non_clique_update_start = Instant::now();
            update_non_clique_samples(conn, &node_ids, &clique_nodes, j, query)?;
            println!("        ✓ Non-clique samples updated [{:.2}ms]",
                     non_clique_update_start.elapsed().as_millis());

            // 7. Delete invalid edges
            let edge_cleanup_start = Instant::now();
            let delete_dist_threshold = if j == levels.len() - 1 {
                1.0
            } else {
                levels[j + 1].1
            };

            delete_edges_by_ids(conn, &clique_edges, delete_dist_threshold)?;
            println!(
                "        🗑️ Deleted {} invalid edges (dist < {:.4}) [{:.2}ms]",
                clique_edges.len(), delete_dist_threshold, edge_cleanup_start.elapsed().as_millis()
            );

            println!("      ✓ Clique level {} completed [{:.2}ms]",
                     j, level_clique_start.elapsed().as_millis());

            let total_clique_time = clique_search_start.elapsed();
            println!("    ✓ Clique search successful [{:.2}ms]", total_clique_time.as_millis());
            return Ok(());
        } else {
            println!("        ○ No clique found at level {} [{:.2}ms]",
                     j, clique_check_start.elapsed().as_millis());
            break;
        }
    }

    let total_clique_time = clique_search_start.elapsed();
    println!("    ○ No cliques from level {} [{:.2}ms]", level, total_clique_time.as_millis());
    Ok(())
}

/// Updates code, code_full and code_state field of `query` at `level`
/// based on data from merged `code` table.
pub fn set_clique_code_for_query(
    conn: &mut Connection,
    query: &mut Query,
    level: usize,
) -> Result<(), duckdb::Error> {
    let level_int_col = format!("L_{}_int", level);

    // 1. Get max code value from `code` table
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

    // 2. Assign code = max_code + 1 and state
    let code_value = max_code_opt.unwrap_or(0) as usize + 1;
    query.code[level] = code_value;
    query.code_state[level] = "C".to_string();

    // 3. Build code_full
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
        "          📝 Query level {} updated • Code: {} • Full: '{}' • State: C",
        level, query.code[level], query.code_full[level]
    );

    Ok(())
}

/// Propagates values from `query` at `level` to all `samples` in the clique
/// in the merged `code` table.
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
        // Use UPSERT (INSERT ... ON CONFLICT DO UPDATE)
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

/// Propagates values from `query` at `level` to nodes NOT in clique
/// with state "S" in merged `code` table.
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
        // Use UPSERT
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

/// Saves all query values to the merged code table
pub fn update_duckdb(
    conn: &mut Connection,
    query: &Query,
) -> DuckResult<()> {
    let sample = &query.sample_name;
    let n = query.code.len();

    // Build dynamic list of columns
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

    // Build dynamic SQL
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

    // Prepare parameters
    let mut params: Vec<Box<dyn duckdb::ToSql>> = Vec::new();
    params.push(Box::new(sample.clone()));

    // Add values in order: int, full, state
    for val in int_values {
        params.push(Box::new(val));
    }
    for val in full_values {
        params.push(Box::new(val.clone()));
    }
    for val in state_values {
        params.push(Box::new(val.clone()));
    }

    // Execute with params_from_iter
    let mut stmt = conn.prepare(&sql)?;
    stmt.execute(duckdb::params_from_iter(params))?;

    Ok(())
}
