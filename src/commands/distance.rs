use clap::Args;
use anyhow::{Result, Context, anyhow};
use hashbrown::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use std::sync::Arc;
use duckdb::Connection;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;

use crate::sketch::sketching::*;

/// Arguments for the distance subcommand
#[derive(Args, Debug)]
pub struct DistanceArgs {
    /// Path to DuckDB database
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,

    /// Path to text file with sample IDs (one per line)
    #[arg(long, required = true, value_name = "IDS_FILE")]
    pub ids: String,

    /// Output file for distance matrix
    #[arg(long, required = true, value_name = "OUTPUT_FILE")]
    pub output: String,

    /// Output format: tsv, csv, or phylip
    #[arg(long, default_value = "tsv", value_name = "FORMAT")]
    pub format: String,

    /// Minimum ANI threshold (0.0-1.0). Only pairs with ANI >= threshold are reported.
    /// For TSV/CSV: outputs only pairs meeting threshold in long format.
    /// For PHYLIP: includes all pairs but useful for filtering before phylogenetic analysis.
    #[arg(long, value_name = "THRESHOLD")]
    pub threshold: Option<f64>,

    /// Number of CPUs to parallelize with rayon
    #[arg(long, value_name = "NCPUS", default_value_t = 1)]
    pub cpus: usize,

    /// Enable verbose output
    #[arg(long)]
    pub verbose: bool,
}

/// Context structure for distance calculation
pub struct DistanceCtx<'a> {
    pub conn: &'a Connection,
    pub metadata_map: HashMap<String, String>,
}

impl<'a> DistanceCtx<'a> {
    pub fn kmer_size(&self) -> usize {
        self.metadata_map["kmer_size"].parse().unwrap_or(21)
    }

    pub fn sketch_size(&self) -> usize {
        self.metadata_map["sketch_size"].parse().unwrap_or(1000)
    }

    pub fn genus(&self) -> &str {
        &self.metadata_map["genus"]
    }

    pub fn acronym(&self) -> &str {
        &self.metadata_map["acronym"]
    }
}

/// Structure to hold a pairwise distance result
#[derive(Debug, Clone)]
pub struct PairwiseDistance {
    pub sample1: String,
    pub sample2: String,
    pub ani: f64,
}

/// Enum to represent where a sample is located
#[derive(Debug, Clone)]
enum SampleLocation {
    InSketches { signature: u64 },
    InDuplicates { original_signature: u64, original_name: String },
}

/// Load metadata from database into HashMap
fn load_metadata_map(conn: &Connection) -> Result<HashMap<String, String>> {
    let start = Instant::now();
    
    // Get column names dynamically
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

/// Normalize a sample name by removing hidden/control characters
fn normalize_sample_name(raw_name: &str) -> String {
    raw_name
        .trim()                       // Elimina espacios, \n, \r, \t
        .trim_matches('\u{FEFF}')     // Elimina BOM (Byte Order Mark)
        .chars()
        .filter(|c| !c.is_control())  // Elimina caracteres de control
        .collect::<String>()
}

/// Read sample names from text file and normalize them
fn read_sample_names(ids_file: &str) -> Result<Vec<String>> {
    let start = Instant::now();
    
    let content = fs::read_to_string(ids_file)
        .with_context(|| format!("Error reading IDs file: {}", ids_file))?;
    
    let names: Vec<String> = content
        .lines()
        .map(|line| normalize_sample_name(line))
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    
    if names.is_empty() {
        return Err(anyhow!("IDs file is empty or contains no valid sample names"));
    }
    
    let duration = start.elapsed();
    println!("✓ Loaded {} sample names from file [{:.2}ms]",
             names.len(), duration.as_millis());
    
    Ok(names)
}

/// Verify sample locations using batch queries (OPTIMIZED)
/// Handles samples in sketches table, duplicates table, or missing
fn map_samples_to_signatures(
    conn: &Connection,
    sample_names: &[String]
) -> Result<HashMap<String, SampleLocation>> {
    let start = Instant::now();
    let mut sample_map = HashMap::new();
    
    // Build parameterized query for batch lookup in sketches
    let placeholders = sample_names.iter()
        .map(|_| "?")
        .collect::<Vec<_>>()
        .join(", ");
    
    let sketches_query = format!(
        "SELECT name, signature FROM sketches WHERE name IN ({})",
        placeholders
    );
    
    // Query all samples from sketches table in one go
    let mut stmt = conn.prepare(&sketches_query)
        .context("Error preparing batch sketches query")?;
    
    let params: Vec<&dyn duckdb::ToSql> = sample_names.iter()
        .map(|s| s as &dyn duckdb::ToSql)
        .collect();
    
    let mut rows = stmt.query(&params[..])
        .context("Error executing batch sketches query")?;
    
    let mut found_in_sketches = HashMap::new();
    while let Some(row) = rows.next()? {
        let name: String = row.get(0)?;
        let sig: u64 = row.get(1)?;
        found_in_sketches.insert(name.clone(), sig);
        sample_map.insert(
            name,
            SampleLocation::InSketches { signature: sig }
        );
    }
    
    // Identify samples not found in sketches
    let not_in_sketches: Vec<&String> = sample_names.iter()
        .filter(|name| !found_in_sketches.contains_key(*name))
        .collect();
    
    if !not_in_sketches.is_empty() {
        // Batch query duplicates table
        let dup_placeholders = not_in_sketches.iter()
            .map(|_| "?")
            .collect::<Vec<_>>()
            .join(", ");
        
        let duplicates_query = format!(
            "SELECT d.sample, d.signature, s.name 
             FROM duplicates d 
             JOIN sketches s ON d.signature = s.signature 
             WHERE d.sample IN ({})",
            dup_placeholders
        );
        
        let mut dup_stmt = conn.prepare(&duplicates_query)
            .context("Error preparing batch duplicates query")?;
        
        let dup_params: Vec<&dyn duckdb::ToSql> = not_in_sketches.iter()
            .map(|s| *s as &dyn duckdb::ToSql)
            .collect();
        
        let mut dup_rows = dup_stmt.query(&dup_params[..])
            .context("Error executing batch duplicates query")?;
        
        let mut found_in_duplicates = HashMap::new();
        while let Some(row) = dup_rows.next()? {
            let sample_name: String = row.get(0)?;
            let original_sig: u64 = row.get(1)?;
            let original_name: String = row.get(2)?;
            
            found_in_duplicates.insert(sample_name.clone(), ());
            sample_map.insert(
                sample_name,
                SampleLocation::InDuplicates {
                    original_signature: original_sig,
                    original_name,
                }
            );
        }
        
        // Check for missing samples
        let missing_names: Vec<String> = not_in_sketches.iter()
            .filter(|name| !found_in_duplicates.contains_key(**name))
            .map(|s| (*s).clone())
            .collect();
        
        if !missing_names.is_empty() {
            for name in &missing_names {
                eprintln!("⚠️  Sample not found in DB: {:?}", name);
            }
            return Err(anyhow!(
                "The following {} sample names were not found in sketches or duplicates tables: {}",
                missing_names.len(),
                missing_names.join(", ")
            ));
        }
    }
    
    let duration = start.elapsed();
    
    // Count samples by location
    let in_sketches = sample_map.values()
        .filter(|loc| matches!(loc, SampleLocation::InSketches { .. }))
        .count();
    let in_duplicates = sample_map.len() - in_sketches;
    
    println!("✓ Verified all {} samples [{:.2}ms]", sample_names.len(), duration.as_millis());
    println!("  - {} in sketches table", in_sketches);
    println!("  - {} in duplicates table", in_duplicates);
    
    Ok(sample_map)
}

/// Extract unique signatures needed for distance calculations
/// Returns a map from signature to the sample names that use it
fn get_unique_signatures_map(
    sample_map: &HashMap<String, SampleLocation>
) -> HashMap<u64, Vec<String>> {
    let mut sig_to_samples: HashMap<u64, Vec<String>> = HashMap::new();
    
    for (sample_name, location) in sample_map {
        let signature = match location {
            SampleLocation::InSketches { signature } => *signature,
            SampleLocation::InDuplicates { original_signature, .. } => *original_signature,
        };
        
        sig_to_samples.entry(signature)
            .or_insert_with(Vec::new)
            .push(sample_name.clone());
    }
    
    sig_to_samples
}

/// Load sketches into memory wrapped in Arc for efficient sharing
fn load_sketches_optimized(
    sketch_manager: &SketchManager,
    signatures: &[u64],
    verbose: bool,
) -> Result<HashMap<u64, Arc<Sketch>>> {
    let start = Instant::now();
    
    let sketch_map: HashMap<u64, Arc<Sketch>> = signatures
        .iter()
        .filter_map(|&sig| {
            sketch_manager.get_sketch(sig)
                .map(|sketch| (sig, Arc::new(sketch.clone())))
        })
        .collect();
    
    let duration = start.elapsed();
    
    if verbose {
        println!("✓ Loaded {} sketches into memory [{:.2}ms]",
                 sketch_map.len(), duration.as_millis());
    }
    
    Ok(sketch_map)
}

/// Calculate pairwise distances in long format (OPTIMIZED)
fn calculate_pairwise_distances_long(
    sketch_map: &HashMap<u64, Arc<Sketch>>,
    sample_names: &[String],
    sample_map: &HashMap<String, SampleLocation>,
    threshold: Option<f64>,
    verbose: bool,
) -> Result<Vec<PairwiseDistance>> {
    let n = sample_names.len();
    
    let distances: Vec<PairwiseDistance> = (0..n)
        .into_par_iter()
        .flat_map(|i| {
            let mut results: Vec<PairwiseDistance> = Vec::new();
            
            let loc_i = sample_map.get(&sample_names[i]).expect("Sample not in map");
            let sig_i = match loc_i {
                SampleLocation::InSketches { signature } => *signature,
                SampleLocation::InDuplicates { original_signature, .. } => *original_signature,
            };
            
            let sketch_i = sketch_map.get(&sig_i).expect("Sketch not found");

            for j in (i + 1)..n {
                let loc_j = sample_map.get(&sample_names[j]).expect("Sample not in map");
                let sig_j = match loc_j {
                    SampleLocation::InSketches { signature } => *signature,
                    SampleLocation::InDuplicates { original_signature, .. } => *original_signature,
                };
                
                let sketch_j = sketch_map.get(&sig_j).expect("Sketch not found");
                let ani = sketch_i.distance(sketch_j);

                if threshold.map_or(true, |thr| ani >= thr) {
                    results.push(PairwiseDistance {
                        sample1: sample_names[i].clone(),
                        sample2: sample_names[j].clone(),
                        ani,
                    });
                }
            }

            if verbose && (i + 1) % 100 == 0 {
                println!("  Progress: {}/{} samples processed", i + 1, n);
            }

            results
        })
        .collect();

    Ok(distances)
}

/// Calculate pairwise distances as full matrix (OPTIMIZED)
fn calculate_pairwise_distances_matrix(
    sketch_map: &HashMap<u64, Arc<Sketch>>,
    sample_names: &[String],
    sample_map: &HashMap<String, SampleLocation>,
    verbose: bool,
) -> Result<Vec<Vec<f64>>> {
    let n = sample_names.len();

    let mut matrix: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut row = vec![0.0f64; n];

            let loc_i = sample_map.get(&sample_names[i]).expect("Sample not in map");
            let sig_i = match loc_i {
                SampleLocation::InSketches { signature } => *signature,
                SampleLocation::InDuplicates { original_signature, .. } => *original_signature,
            };
            
            let sketch_i = sketch_map.get(&sig_i).expect("Sketch not found");

            for j in 0..n {
                if i == j {
                    row[j] = 1.0;
                } else if i < j {
                    let loc_j = sample_map.get(&sample_names[j]).expect("Sample not in map");
                    let sig_j = match loc_j {
                        SampleLocation::InSketches { signature } => *signature,
                        SampleLocation::InDuplicates { original_signature, .. } => *original_signature,
                    };
                    
                    let sketch_j = sketch_map.get(&sig_j).expect("Sketch not found");
                    row[j] = sketch_i.distance(sketch_j);
                }
            }

            if verbose && (i + 1) % 100 == 0 {
                println!("  Progress: {}/{} samples processed", i + 1, n);
            }

            row
        })
        .collect();

    // Symmetrize lower triangle
    for i in 0..n {
        for j in 0..i {
            matrix[i][j] = matrix[j][i];
        }
    }

    Ok(matrix)
}

/// Save distances in TSV long format
fn save_distances_tsv_long(
    output_file: &str,
    distances: &[PairwiseDistance],
) -> Result<()> {
    let start = Instant::now();
    let file = File::create(output_file)
        .with_context(|| format!("Error creating output file: {}", output_file))?;
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "Sample1\tSample2\tANI\tANI_percent\tDistance")?;
    
    for dist in distances {
        writeln!(
            writer,
            "{}\t{}\t{:.6}\t{:.2}\t{:.6}",
            dist.sample1,
            dist.sample2,
            dist.ani,
            dist.ani * 100.0,
            1.0 - dist.ani
        )?;
    }
    
    writer.flush()?;
    
    let duration = start.elapsed();
    println!("✓ Distances saved as TSV (long format) [{:.2}ms]", duration.as_millis());
    
    Ok(())
}

/// Save distances in CSV long format
fn save_distances_csv_long(
    output_file: &str,
    distances: &[PairwiseDistance],
) -> Result<()> {
    let start = Instant::now();
    let file = File::create(output_file)
        .with_context(|| format!("Error creating output file: {}", output_file))?;
    let mut writer = BufWriter::new(file);
    
    writeln!(writer, "Sample1,Sample2,ANI,ANI_percent,Distance")?;
    
    for dist in distances {
        writeln!(
            writer,
            "{},{},{:.6},{:.2},{:.6}",
            dist.sample1,
            dist.sample2,
            dist.ani,
            dist.ani * 100.0,
            1.0 - dist.ani
        )?;
    }
    
    writer.flush()?;
    
    let duration = start.elapsed();
    println!("✓ Distances saved as CSV (long format) [{:.2}ms]", duration.as_millis());
    
    Ok(())
}

/// Save distance matrix in PHYLIP format
fn save_matrix_phylip(
    output_file: &str,
    sample_names: &[String],
    matrix: &[Vec<f64>],
) -> Result<()> {
    let start = Instant::now();
    let file = File::create(output_file)
        .with_context(|| format!("Error creating output file: {}", output_file))?;
    let mut writer = BufWriter::new(file);
    
    let n = sample_names.len();
    
    writeln!(writer, "{}", n)?;
    
    for (i, row) in matrix.iter().enumerate() {
        let name = if sample_names[i].len() > 10 {
            &sample_names[i][..10]
        } else {
            &sample_names[i]
        };
        write!(writer, "{:<10}", name)?;
        
        for &value in row {
            let distance = 1.0 - value;
            write!(writer, " {:.6}", distance)?;
        }
        writeln!(writer)?;
    }
    
    writer.flush()?;
    
    let duration = start.elapsed();
    println!("✓ Distance matrix saved as PHYLIP [{:.2}ms]", duration.as_millis());
    
    Ok(())
}

/// Calculate statistics from distances in long format
fn calculate_statistics_long(distances: &[PairwiseDistance]) -> Option<DistanceStats> {
    if distances.is_empty() {
        return None;
    }
    
    let mut anis: Vec<f64> = distances.iter().map(|d| d.ani).collect();
    anis.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let min_ani = anis[0];
    let max_ani = anis[anis.len() - 1];
    let mean_ani = anis.iter().sum::<f64>() / anis.len() as f64;
    let median_ani = anis[anis.len() / 2];
    
    Some(DistanceStats {
        count: distances.len(),
        min_ani,
        max_ani,
        mean_ani,
        median_ani,
    })
}

/// Calculate statistics from distance matrix
fn calculate_statistics_matrix(matrix: &[Vec<f64>]) -> Option<DistanceStats> {
    let mut all_distances = Vec::new();
    
    for i in 0..matrix.len() {
        for j in (i+1)..matrix[i].len() {
            all_distances.push(matrix[i][j]);
        }
    }
    
    if all_distances.is_empty() {
        return None;
    }
    
    all_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    
    let min_ani = all_distances[0];
    let max_ani = all_distances[all_distances.len() - 1];
    let mean_ani = all_distances.iter().sum::<f64>() / all_distances.len() as f64;
    let median_ani = all_distances[all_distances.len() / 2];
    
    Some(DistanceStats {
        count: all_distances.len(),
        min_ani,
        max_ani,
        mean_ani,
        median_ani,
    })
}

/// Structure to hold distance statistics
#[derive(Debug)]
struct DistanceStats {
    count: usize,
    min_ani: f64,
    max_ani: f64,
    mean_ani: f64,
    median_ani: f64,
}

impl DistanceStats {
    fn print(&self) {
        println!("\n=== Distance Statistics ===");
        println!("Number of pairwise comparisons: {}", self.count);
        println!("Minimum ANI: {:.6} ({:.2}%)", self.min_ani, self.min_ani * 100.0);
        println!("Maximum ANI: {:.6} ({:.2}%)", self.max_ani, self.max_ani * 100.0);
        println!("Mean ANI: {:.6} ({:.2}%)", self.mean_ani, self.mean_ani * 100.0);
        println!("Median ANI: {:.6} ({:.2}%)", self.median_ani, self.median_ani * 100.0);
    }
}

/// Main function for the distance command
pub fn distance_command(args: &DistanceArgs) -> Result<()> {
    let command_start = Instant::now();
    
    println!("=== Starting distance command ===");
    println!("Database: {}", args.db);
    println!("IDs file: {}", args.ids);
    println!("Output file: {}", args.output);
    println!("Format: {}", args.format);
    println!("CPUs: {}", args.cpus);
    
    if let Some(thr) = args.threshold {
        println!("ANI Threshold: >= {:.4} ({:.2}%)", thr, thr * 100.0);
    } else {
        println!("ANI Threshold: None (all pairs will be reported)");
    }

    // Configure Rayon thread pool
    ThreadPoolBuilder::new()
        .num_threads(args.cpus)
        .build_global()
        .context("Error configuring Rayon thread pool")?;
    
    let actual_threads = rayon::current_num_threads();
    println!("✓ Rayon thread pool configured with {} threads", actual_threads);
    
    // Validate input files
    let validation_start = Instant::now();
    if !Path::new(&args.db).exists() {
        return Err(anyhow!("Database not found: {}", args.db));
    }
    if !Path::new(&args.ids).exists() {
        return Err(anyhow!("IDs file not found: {}", args.ids));
    }
    println!("✓ Input files validated [{:.2}ms]", 
             validation_start.elapsed().as_millis());
    
    // Validate output format
    let valid_formats = ["tsv", "csv", "phylip"];
    if !valid_formats.contains(&args.format.as_str()) {
        return Err(anyhow!(
            "Invalid format '{}'. Valid formats: {}",
            args.format,
            valid_formats.join(", ")
        ));
    }
    
    // Validate threshold if provided
    if let Some(thr) = args.threshold {
        if thr < 0.0 || thr > 1.0 {
            return Err(anyhow!(
                "Invalid threshold '{}'. Must be between 0.0 and 1.0",
                thr
            ));
        }
        
        if args.format == "phylip" {
            println!("⚠️  Warning: Threshold is ignored for PHYLIP format (full matrix required)");
        }
    }
    
    // Open database connection
    let db_connect_start = Instant::now();
    let conn = Connection::open(&args.db)
        .with_context(|| format!("Error opening database: {}", args.db))?;
    println!("✓ Database connection established [{:.2}ms]",
             db_connect_start.elapsed().as_millis());
    
    // Load metadata
    let metadata_map = load_metadata_map(&conn)?;
    
    let ctx = DistanceCtx {
        conn: &conn,
        metadata_map,
    };
    
    println!("\n=== Project information ===");
    println!("Genus: {}", ctx.genus());
    println!("Acronym: {}", ctx.acronym());
    println!("K-mer size: {}", ctx.kmer_size());
    println!("Sketch size: {}", ctx.sketch_size());
    
    // Read sample names from file
    let sample_names = read_sample_names(&args.ids)?;
    
    // Map samples to their signatures (OPTIMIZED: batch queries)
    let sample_map = map_samples_to_signatures(&conn, &sample_names)?;
    
    // Load SketchManager from database
    let sketch_load_start = Instant::now();
    let sketch_manager = load_sketch_manager_from_db(
        &conn,
        ctx.kmer_size(),
        ctx.sketch_size()
    ).context("Error loading SketchManager from database")?;
    
    println!("✓ SketchManager loaded from DB. Contains {} sketches [{:.2}s]",
             sketch_manager.length(),
             sketch_load_start.elapsed().as_secs_f32());
    
    // Get unique signatures and load sketches into memory
    let sig_to_samples = get_unique_signatures_map(&sample_map);
    let unique_signatures: Vec<u64> = sig_to_samples.keys().copied().collect();
    
    let sketch_map = load_sketches_optimized(
        &sketch_manager,
        &unique_signatures,
        args.verbose
    )?;
    
    // Verify all required signatures are loaded
    let missing_sigs: Vec<u64> = unique_signatures
        .iter()
        .copied()
        .filter(|sig| !sketch_map.contains_key(sig))  
        .collect();

    if !missing_sigs.is_empty() {
        return Err(anyhow!(
            "The following {} signatures have no sketches: {:?}",
            missing_sigs.len(),
            missing_sigs
        ));
    }
    
    // Calculate distances and save based on format
    match args.format.as_str() {
        "tsv" => {
            let distances = calculate_pairwise_distances_long(
                &sketch_map,
                &sample_names,
                &sample_map,
                args.threshold,
                args.verbose
            )?;
            
            if let Some(stats) = calculate_statistics_long(&distances) {
                stats.print();
            }
            
            save_distances_tsv_long(&args.output, &distances)?;
        }
        
        "csv" => {
            let distances = calculate_pairwise_distances_long(
                &sketch_map,
                &sample_names,
                &sample_map,
                args.threshold,
                args.verbose
            )?;
            
            if let Some(stats) = calculate_statistics_long(&distances) {
                stats.print();
            }
            
            save_distances_csv_long(&args.output, &distances)?;
        }
        
        "phylip" => {
            let matrix = calculate_pairwise_distances_matrix(
                &sketch_map,
                &sample_names,
                &sample_map,
                args.verbose
            )?;
            
            if let Some(stats) = calculate_statistics_matrix(&matrix) {
                stats.print();
            }
            
            save_matrix_phylip(&args.output, &sample_names, &matrix)?;
        }
        
        _ => unreachable!(),
    }
    
    let total_duration = command_start.elapsed();
    println!("\n=== Distance command completed successfully [{:.2}s] ===",
             total_duration.as_secs_f32());
    
    Ok(())
}
