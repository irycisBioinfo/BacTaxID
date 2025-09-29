
use clap::Args;
use anyhow::{Result, Context, anyhow};
use hashbrown::HashMap;
use std::fs;
use std::path::Path;
use std::time::Instant;
use std::io::Write;
use duckdb::{Connection, Row, ToSql, params, Result as DuckResult};
use crate::{
    sketch::sketching::*,
    db::db::*,
};

/// Arguments for the classify subcommand
#[derive(Args, Debug)]
pub struct ClassifyArgs {
    /// Path to DuckDB database
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,
    /// Plain text file with paths to FASTA files (one per line)
    #[arg(long, required = true, value_name = "QUERY_FILES")]
    pub queries: String,
    /// Output file for classification results (optional)
    #[arg(long, required = true, value_name = "OUTPUT_FILE")]
    pub output: Option<String>,
    /// Enable verbose output
    #[arg(long)]
    pub verbose: bool,
}

/// Classification result for a single query
#[derive(Debug, Clone)]
pub struct ClassificationResult {
    pub query_id: String,
    pub best_hit_id: String,
    pub similarity_score: f64,
    pub levels_reached: usize,
    pub final_code: Vec<usize>,
    pub final_code_full: Vec<String>,
    pub best_hit_code: Vec<usize>,
    pub best_hit_code_full: Vec<String>,
}

/// Context for classification (reusing UpdateCtx structure)
pub struct ClassifyCtx<'a> {
    /// Mutable connection to DuckDB (for queries only, no insertions)
    pub conn: &'a mut Connection,
    /// Levels table (L_0 ... L_N → dist) as ordered vector
    pub levels: Vec<(String, f64)>,
    /// Entire metadata row as (column → value in UTF-8)  
    pub metadata_map: HashMap<String, String>,
}

impl<'a> ClassifyCtx<'a> {
    /// Convenience accessors for metadata (same as UpdateCtx)
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
    
    pub fn reference_size(&self) -> usize {
        self.metadata_map["reference_size"].parse().unwrap_or(100)
    }
    
    /// Number of levels
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }
    
    /// Access to mutable connection
    pub fn connection_mut(&mut self) -> &mut Connection {
        self.conn
    }
}

/// Query structure for classification (simplified version of update Query)
pub struct ClassifyQuery {
    pub sample_name: String,
    pub code: Vec<usize>,
    pub code_full: Vec<String>,
    pub sketch: SketchManager,
}

impl ClassifyQuery {
    /// Create a new ClassifyQuery (similar to update Query::new but without code_state)
    pub fn new(path: &Path, ctx: &ClassifyCtx) -> Result<Self> {
        let query_creation_start = Instant::now();
        let sample_name = path.file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();
        let n_levels = ctx.levels.len();
        
        // Initialize vectors
        let code = vec![0_usize; n_levels];
        let code_full = vec!["".to_string(); n_levels];
        
        // Create sketch
        let sketch_creation_start = Instant::now();
        let mut sketch_manager = SketchManager::new(ctx.kmer_size(), ctx.sketch_size());
        let sketch = Sketch::new(path, ctx.kmer_size(), ctx.sketch_size())?;
        sketch_manager.add_sketch(sketch)?;
        
        if ctx.metadata_map.get("verbose").unwrap_or(&"false".to_string()) == "true" {
            println!(" ✓ Sketch created [{:.2}ms]", sketch_creation_start.elapsed().as_millis());
            println!(" ✓ Query initialized [{:.2}ms]", query_creation_start.elapsed().as_millis());
        }
        
        Ok(ClassifyQuery {
            sample_name,
            code,
            code_full,
            sketch: sketch_manager,
        })
    }
}

/// Load levels table (exact copy from update.rs)
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

/// Load metadata map (exact copy from update.rs)
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

/// Retrieve classifiers from the merged code table (exact copy from update.rs)
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
        // Unrecognized condition
        _ => {
            return Err(anyhow!("Invalid condition '{}' for classify. Only 'C' is supported", condition));
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
        println!(" └─ SQL query executed [{:.2}ms]", duration.as_millis());
    }
    
    Ok(result)
}

/// Find best hit (exact copy from update.rs)
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

/// Classify a single query file (simplified version of update_single_file)
/// DETIENE LA BÚSQUEDA cuando no se encuentra un best hit sobre los classifiers
pub fn classify_single_file(
    fasta_path: &Path,
    ctx: &mut ClassifyCtx,
    sketch_manager: &SketchManager,
    verbose: bool
) -> Result<ClassificationResult> {
    let file_processing_start = Instant::now();
    
    // 1. Create Query (includes validations, sketch and vectors)
    let query_start = Instant::now();
    let mut query = ClassifyQuery::new(fasta_path, ctx)?;
    
    if verbose {
        println!(" ✓ Query created • sample: {} • k: {} • sketch: {} • levels: {} [{:.2}ms]",
            query.sample_name,
            ctx.kmer_size(),
            ctx.sketch_size(),
            ctx.num_levels(),
            query_start.elapsed().as_millis()
        );
    }
    
    // 2. Processing variables
    let mut ref_db: Vec<(String, usize, String, String)>;
    let mut distances: Vec<(String, String, f64)>;
    let mut bh: Option<(String, usize, String, String)> = None;
    let mut ref_ids: Vec<String> = Vec::new();
    let mut levels_reached = 0;
    let mut final_best_hit = None;
    let mut final_similarity = 0.0;
    
    // 3. Process each level in ascending order of distance
    for i in 0..ctx.levels.len() {
        let level_start = Instant::now();
        
        if verbose {
            println!(" 🔍 Processing level {} [{:.4} threshold]", i, ctx.levels[i].1);
        }
        
        // 3.1 Retrieve classifiers with condition "C" ONLY
        let classifiers_start = Instant::now();
        ref_db = if i == 0 {
            retrieve_classifiers(ctx.connection_mut(), i, "", "C")?
        } else {
            retrieve_classifiers(ctx.connection_mut(), i, &query.code_full[i - 1], "C")?
        };
        ref_ids = ref_db.iter().map(|(sample, _, _, _)| sample.clone()).collect();
        
        if verbose {
            println!(" ✓ Classifiers obtained: {} [{:.2}ms]",
                ref_ids.len(), classifiers_start.elapsed().as_millis());
        }
        
        // 3.2 Si no hay classifiers, detener la búsqueda
        if ref_db.is_empty() {
            if verbose {
                println!(" ○ No classifiers found at level {}, stopping classification [{:.2}ms]", 
                    i, level_start.elapsed().as_millis());
            }
            break;
        }
        
        // 3.3 Calculate distances
        let distances_start = Instant::now();
        distances = pw_one_to_many(
            &query.sketch,
            sketch_manager,
            &ref_ids,
            ctx.levels[i].1,
        );
        
        if verbose {
            println!(" ✓ Distances calculated: {} [{:.2}ms]",
                distances.len(), distances_start.elapsed().as_millis());
        }
        
        // 3.4 Si no hay distancias válidas, detener la búsqueda
        if distances.is_empty() {
            if verbose {
                println!(" ○ No valid distances at level {}, stopping classification [{:.2}ms]", 
                    i, level_start.elapsed().as_millis());
            }
            break;
        }
        
        // 3.5 Find best hit
        let best_hit_start = Instant::now();
        bh = best_hit(&distances, &ref_db);
        
        if verbose {
            println!(" ✓ Best hit calculated [{:.2}ms]", best_hit_start.elapsed().as_millis());
        }
        
        // 3.6 Si no hay best hit, detener la búsqueda
        let Some(best_hit_tuple) = bh else {
            if verbose {
                println!(" ○ No best hit found at level {}, stopping classification [{:.2}ms]", 
                    i, level_start.elapsed().as_millis());
            }
            break;
        };
        
        // 3.7 Process best hit
        let code_val = best_hit_tuple.1;
        let code_full_val = best_hit_tuple.2.clone();
        
        // 3.8 Si el código es 0, detener
        if code_val == 0 {
            if verbose {
                println!(" ○ Best hit has code = 0 at level {}, stopping classification [{:.2}ms]", 
                    i, level_start.elapsed().as_millis());
            }
            break;
        }
        
        // 3.9 Update query with valid best hit
        query.code[i] = code_val;
        query.code_full[i] = code_full_val.clone();
        
        // 3.10 Consistency check with previous level
        if i > 0 && !query.code_full[i].starts_with(&query.code_full[i-1]) {
            if verbose {
                println!(" ○ Code_full inconsistency at level {} ('{}' does not start with '{}'), stopping classification [{:.2}ms]", 
                    i, query.code_full[i], query.code_full[i-1], level_start.elapsed().as_millis());
            }
            break;
        }
        
        // 3.11 Calculate similarity for this best hit
        if let Some(best_distance_tuple) = distances.iter().find(|(_, ref_name, _)| ref_name == &best_hit_tuple.0) {
            final_similarity = best_distance_tuple.2;
        }
        
        // 3.12 Store current best hit as final result
        final_best_hit = Some(best_hit_tuple);
        levels_reached = i + 1;
        
        if verbose {
            println!(" ✓ Level {} completed successfully • code: {} • full: '{}' • sim: {:.4} [{:.2}ms]", 
                i, code_val, code_full_val, final_similarity, level_start.elapsed().as_millis());
        }
        
        // Continue to next level (no early stopping conditions in classify)
    }
    
    // 4. Create classification result
    let best_hit_info = final_best_hit.unwrap_or(("Unknown".to_string(), 0, "".to_string(), "".to_string()));
    
    // Get best hit full classification from database
    let mut best_hit_code = vec![0_usize; ctx.num_levels()];
    let mut best_hit_code_full = vec!["".to_string(); ctx.num_levels()];
    
    // Query for best hit classification
    let num_levels = ctx.num_levels();
    if let Ok(hit_classification) = get_sample_full_classification(ctx.connection_mut(), &best_hit_info.0, num_levels) {
        best_hit_code = hit_classification.0;
        best_hit_code_full = hit_classification.1;
    }
    
    let result = ClassificationResult {
        query_id: query.sample_name.clone(),
        best_hit_id: best_hit_info.0,
        similarity_score: final_similarity,
        levels_reached,
        final_code: query.code,
        final_code_full: query.code_full,
        best_hit_code,
        best_hit_code_full,
    };
    
    if verbose {
        println!(
            "✓ Classification complete • sample: {} • levels_reached: {} • best_hit: {} • similarity: {:.4} [{:.2}s]",
            result.query_id, result.levels_reached, result.best_hit_id, result.similarity_score, 
            file_processing_start.elapsed().as_secs_f32()
        );
    }
    
    Ok(result)
}

/// Get full classification for a sample from database
fn get_sample_full_classification(
    conn: &mut Connection, 
    sample_name: &str, 
    num_levels: usize
) -> Result<(Vec<usize>, Vec<String>)> {
    let mut int_cols = Vec::new();
    let mut full_cols = Vec::new();
    
    for i in 0..num_levels {
        int_cols.push(format!("L_{}_int", i));
        full_cols.push(format!("L_{}_full", i));
    }
    
    let all_cols = [&int_cols[..], &full_cols[..]].concat();
    let sql = format!("SELECT {} FROM code WHERE sample = ? LIMIT 1", all_cols.join(", "));
    
    let mut stmt = conn.prepare(&sql)?;
    let row = stmt.query_row([sample_name], |row| {
        let mut codes = Vec::new();
        let mut code_fulls = Vec::new();
        
        // Read int values
        for i in 0..num_levels {
            let code: Option<i64> = row.get(i)?;
            codes.push(code.unwrap_or(0) as usize);
        }
        
        // Read full values
        for i in 0..num_levels {
            let code_full: Option<String> = row.get(num_levels + i)?;
            code_fulls.push(code_full.unwrap_or_default());
        }
        
        Ok((codes, code_fulls))
    })?;
    
    Ok(row)
}



fn save_results(results: &[ClassificationResult], output_path: &str) -> anyhow::Result<()> {
    let mut output_file = std::fs::File::create(output_path)?;
    // Cabecera
    writeln!(output_file, "query_id\tbest_hit_id\tsimilarity_score\tlevels_reached\tfinal_code\tbest_hit_code")?;
    for result in results {
        writeln!(
            output_file,
            "{}\t{}\t{:.6}\t{}\t{}\t{}",
            result.query_id,
            result.best_hit_id,
            result.similarity_score,
            result.levels_reached,
            result.final_code.iter().map(|v| v.to_string()).collect::<Vec<_>>().join("."),
            result.best_hit_code.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(".")            
        )?;
    }
    Ok(())
}
/// Main function for the classify command
pub fn classify_command(args: &ClassifyArgs) -> Result<()> {
    let command_start = Instant::now();
    println!("=== Starting classify command (read-only mode) ===");
    println!("Database: {}", args.db);
    println!("Queries file: {}", args.queries);
    
    if let Some(ref output) = args.output {
        println!("Output file: {}", output);
    }
    
    // Check that files exist
    let validation_start = Instant::now();
    if !Path::new(&args.db).exists() {
        return Err(anyhow!("Database not found: {}", args.db));
    }
    if !Path::new(&args.queries).exists() {
        return Err(anyhow!("Queries file not found: {}", args.queries));
    }
    println!("✓ Files validated [{:.2}ms]", validation_start.elapsed().as_millis());
    
    // Open mutable connection to database (read-only)
    let db_connect_start = Instant::now();
    let mut conn = Connection::open(&args.db)
        .with_context(|| format!("Error opening database: {}", args.db))?;
    println!("✓ Database connection established [{:.2}ms]",
        db_connect_start.elapsed().as_millis());
    
    // Preload levels and metadata tables
    let levels = load_levels_vec(&conn)?;
    let mut metadata_map = load_metadata_map(&conn)?;
    
    // Add verbose flag to metadata for internal use
    metadata_map.insert("verbose".to_string(), args.verbose.to_string());
    
    // Create shared context
    let mut ctx = ClassifyCtx {
        conn: &mut conn,
        levels,
        metadata_map,
    };
    
    println!("\n=== Project information ===");
    println!("Genus: {}", ctx.genus());
    println!("Acronym: {}", ctx.acronym());
    println!("K-mer size: {}", ctx.kmer_size());
    println!("Sketch size: {}", ctx.sketch_size());
    println!("Available levels: {}", ctx.num_levels());
    
    // Load SketchManager from DuckDB
    let sketch_load_start = Instant::now();
    let sketch_manager_result = load_sketch_manager_from_db(
        ctx.conn,
        ctx.kmer_size(),
        ctx.sketch_size()
    );
    
    println!("✓ SketchManager loaded from DB [{:.2}ms]",
        sketch_load_start.elapsed().as_millis());
    
    match sketch_manager_result {
        Ok(sketch_manager) => {
            println!(
                "✓ SketchManager loaded from DB. Contains {} sketches",
                sketch_manager.length()
            );
            
            // Read file list
            let file_read_start = Instant::now();
            let file_content = fs::read_to_string(&args.queries)
                .with_context(|| format!("Error reading queries file: {}", args.queries))?;
            
            let files: Vec<&str> = file_content
                .lines()
                .map(|line| line.trim())
                .filter(|line| !line.is_empty() && !line.starts_with('#'))
                .collect();
            
            println!("✓ Query file list read: {} files [{:.2}ms]",
                files.len(), file_read_start.elapsed().as_millis());
            
            // Process each file
            let mut results = Vec::new();
            let mut processed_files = 0;
            let mut skipped_files = 0;
            let processing_start = Instant::now();
            
            for line in files {
                let fasta_path = Path::new(line);
                if !fasta_path.exists() {
                    eprintln!("FASTA file not found: {}", fasta_path.display());
                    skipped_files += 1;  
                    continue;
                }
                
                let file_start = Instant::now();
                match classify_single_file(fasta_path, &mut ctx, &sketch_manager, args.verbose) {
                    Ok(result) => {
                        results.push(result);
                        processed_files += 1;
                        
                        if args.verbose {
                            println!(" └─ File processed [{:.2}s]", file_start.elapsed().as_secs_f32());
                        }
                    }
                    Err(e) => {
                        eprintln!("Error processing {}: {}", fasta_path.display(), e);
                        skipped_files += 1;
                    }
                }
            }
            
            println!("✓ File processing completed: {} processed, {} skipped [{:.2}s]",
                processed_files, skipped_files, processing_start.elapsed().as_secs_f32());
            
            // Save results if output file specified
            if let Some(ref output_path) = args.output {
                save_results(&results, output_path)?;
            }
            
            // Print summary
            println!("\n=== Classification Summary ===");
            println!("Total queries processed: {}", results.len());
            
            if args.verbose && !results.is_empty() {
                let avg_levels = results.iter().map(|r| r.levels_reached).sum::<usize>() as f64 / results.len() as f64;
                let avg_similarity = results.iter().map(|r| r.similarity_score).sum::<f64>() / results.len() as f64;
                
                println!("Average levels reached: {:.2}", avg_levels);
                println!("Average similarity: {:.4}", avg_similarity);
                
                // Distribution by levels reached
                let mut level_counts = HashMap::new();
                for result in &results {
                    *level_counts.entry(result.levels_reached).or_insert(0) += 1;
                }
                
                println!("\nDistribution by levels reached:");
                let mut sorted_levels: Vec<_> = level_counts.iter().collect();
                sorted_levels.sort_by_key(|&(level, _)| level);
                for (level, count) in sorted_levels {
                    println!("  {} levels: {} queries", level, count);
                }
                
                // Show some examples
                println!("\nExample results:");
                for (i, result) in results.iter().take(5).enumerate() {
                    println!("  {}: {} -> {} (sim: {:.4}, levels: {})", 
                        i + 1,
                        result.query_id, 
                        result.best_hit_id,
                        result.similarity_score,
                        result.levels_reached
                    );
                }
            }
        }
        Err(e) => {
            return Err(anyhow!("Error loading SketchManager from database: {}", e));
        }
    }
    
    let total_duration = command_start.elapsed();
    println!("\n=== Classify command completed successfully [{:.2}s] ===",
        total_duration.as_secs_f32());
    
    Ok(())
}
