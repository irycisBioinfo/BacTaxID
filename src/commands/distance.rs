use clap::Args;
use anyhow::{Result, Context, anyhow};
use hashbrown::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use duckdb::Connection;
use rayon::prelude::*;

use crate::{
    sketch::sketching::*,
    db::db::*,
};

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

    /// Number of threads for parallel computation
    #[arg(long, default_value_t = 1, value_name = "THREADS")]
    pub threads: usize,

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

/// Read sample IDs from text file
fn read_sample_ids(ids_file: &str) -> Result<Vec<String>> {
    let start = Instant::now();
    
    let content = fs::read_to_string(ids_file)
        .with_context(|| format!("Error reading IDs file: {}", ids_file))?;
    
    let ids: Vec<String> = content
        .lines()
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    
    if ids.is_empty() {
        return Err(anyhow!("IDs file is empty or contains no valid IDs"));
    }
    
    let duration = start.elapsed();
    println!("✓ Loaded {} sample IDs from file [{:.2}ms]",
             ids.len(), duration.as_millis());
    
    Ok(ids)
}

/// Verify that all sample IDs exist in the database
fn verify_samples_exist(conn: &Connection, ids: &[String]) -> Result<()> {
    let start = Instant::now();
    let mut missing_ids = Vec::new();
    
    for id in ids {
        let exists: bool = conn.query_row(
            "SELECT COUNT(*) > 0 FROM sketches WHERE sample = ?",
            [id],
            |row| row.get(0)
        )?;
        
        if !exists {
            missing_ids.push(id.clone());
        }
    }
    
    if !missing_ids.is_empty() {
        return Err(anyhow!(
            "The following {} sample IDs were not found in the database: {}",
            missing_ids.len(),
            missing_ids.join(", ")
        ));
    }
    
    let duration = start.elapsed();
    println!("✓ Verified all {} samples exist in database [{:.2}ms]",
             ids.len(), duration.as_millis());
    
    Ok(())
}

/// Calculate pairwise distances in long format with optional threshold
fn calculate_pairwise_distances_long(
    sketch_manager: &SketchManager,
    ids: &[String],
    threshold: Option<f64>,
    threads: usize,
    verbose: bool,
) -> anyhow::Result<Vec<PairwiseDistance>> {
    let n = ids.len();

    // Configurar pool con num_threads y anotar tipo de retorno explícito
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("Error configuring thread pool")?;

    let distances: Vec<PairwiseDistance> = pool.install(|| {
        // Anotar tipo en collect para evitar E0283
        (0..n)
            .into_par_iter()
            .flat_map(|i| {
                // También podemos anotar el vector intermedio
                let mut results: Vec<PairwiseDistance> = Vec::new();
                let sketch_i = sketch_manager
                    .get_sketch(&ids[i])
                    .expect("Sketch not found");

                for j in (i + 1)..n {
                    let sketch_j = sketch_manager
                        .get_sketch(&ids[j])
                        .expect("Sketch not found");
                    let ani = sketch_i.distance(sketch_j);

                    if threshold.map_or(true, |thr| ani >= thr) {
                        results.push(PairwiseDistance {
                            sample1: ids[i].clone(),
                            sample2: ids[j].clone(),
                            ani,
                        });
                    }
                }

                if verbose && (i + 1) % 10 == 0 {
                    println!("  Progress: {}/{} samples processed", i + 1, n);
                }

                results
            })
            .collect::<Vec<PairwiseDistance>>() // <- clave: anotar tipo aquí
    });

    Ok(distances)
}


/// Calculate pairwise distances as full matrix (for PHYLIP format)
fn calculate_pairwise_distances_matrix(
    sketch_manager: &SketchManager,
    ids: &[String],
    threads: usize,
    verbose: bool,
) -> anyhow::Result<Vec<Vec<f64>>> {
    let n = ids.len();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("Error configuring thread pool")?;

    // Tipo explícito del valor devuelto por install
    let mut matrix: Vec<Vec<f64>> = pool.install(|| {
        (0..n)
            .into_par_iter()
            .map(|i| {
                let mut row = vec![0.0f64; n];

                // Obtener referencias inmutables locales (no capturar nada no Send)
                let sketch_i = sketch_manager
                    .get_sketch(&ids[i])
                    .expect("Sketch not found");

                for j in 0..n {
                    if i == j {
                        row[j] = 1.0;
                    } else if i < j {
                        let sketch_j = sketch_manager
                            .get_sketch(&ids[j])
                            .expect("Sketch not found");
                        row[j] = sketch_i.distance(sketch_j);
                    }
                }

                if verbose && (i + 1) % 10 == 0 {
                    println!("  Progress: {}/{} samples processed", i + 1, n);
                }

                row
            })
            .collect::<Vec<Vec<f64>>>() // <- anotar tipo de collect
    });

    // Simetrizar el triángulo inferior
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
    
    // Write header
    writeln!(writer, "Sample1\tSample2\tANI\tANI_percent\tDistance")?;
    
    // Write distances
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
    
    // Write header
    writeln!(writer, "Sample1,Sample2,ANI,ANI_percent,Distance")?;
    
    // Write distances
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
    ids: &[String],
    matrix: &[Vec<f64>],
) -> Result<()> {
    let start = Instant::now();
    let file = File::create(output_file)
        .with_context(|| format!("Error creating output file: {}", output_file))?;
    let mut writer = BufWriter::new(file);
    
    let n = ids.len();
    
    // Write header (number of taxa)
    writeln!(writer, "{}", n)?;
    
    // Write matrix rows
    // PHYLIP format: first 10 characters for name, then distances
    for (i, row) in matrix.iter().enumerate() {
        let name = if ids[i].len() > 10 {
            &ids[i][..10]
        } else {
            &ids[i]
        };
        write!(writer, "{:<10}", name)?;
        
        for &value in row {
            // Convert ANI to distance (1 - ANI)
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
    println!("Threads: {}", args.threads);
    
    if let Some(thr) = args.threshold {
        println!("ANI Threshold: >= {:.4} ({:.2}%)", thr, thr * 100.0);
    } else {
        println!("ANI Threshold: None (all pairs will be reported)");
    }
    
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
        
        // Warn if threshold is used with PHYLIP format
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
    
    // Read sample IDs
    let ids = read_sample_ids(&args.ids)?;
    
    // Verify all samples exist in database
    verify_samples_exist(&conn, &ids)?;
    
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
    
    // Verify all requested samples are in SketchManager
    let missing_sketches: Vec<&String> = ids.iter()
        .filter(|id| !sketch_manager.contains(id))
        .collect();
    
    if !missing_sketches.is_empty() {
        return Err(anyhow!(
            "The following {} samples have no sketches in the database: {}",
            missing_sketches.len(),
            missing_sketches.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
        ));
    }
    
    // Calculate distances and save based on format
    match args.format.as_str() {
        "tsv" => {
            let distances = calculate_pairwise_distances_long(
                &sketch_manager,
                &ids,
                args.threshold,
                args.threads,
                args.verbose
            )?;
            
            if let Some(stats) = calculate_statistics_long(&distances) {
                stats.print();
            }
            
            save_distances_tsv_long(&args.output, &distances)?;
        }
        
        "csv" => {
            let distances = calculate_pairwise_distances_long(
                &sketch_manager,
                &ids,
                args.threshold,
                args.threads,
                args.verbose
            )?;
            
            if let Some(stats) = calculate_statistics_long(&distances) {
                stats.print();
            }
            
            save_distances_csv_long(&args.output, &distances)?;
        }
        
        "phylip" => {
            // PHYLIP requires full matrix, threshold is not applied
            let matrix = calculate_pairwise_distances_matrix(
                &sketch_manager,
                &ids,
                args.threads,
                args.verbose
            )?;
            
            if let Some(stats) = calculate_statistics_matrix(&matrix) {
                stats.print();
            }
            
            save_matrix_phylip(&args.output, &ids, &matrix)?;
        }
        
        _ => unreachable!(), // Already validated
    }
    
    let total_duration = command_start.elapsed();
    println!("\n=== Distance command completed successfully [{:.2}s] ===",
             total_duration.as_secs_f32());
    
    Ok(())
}
