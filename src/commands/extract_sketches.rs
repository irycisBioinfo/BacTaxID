use clap::Args;
use anyhow::{Result, Context};
use std::fs::File;
use std::path::Path;
use std::io::BufWriter;
use duckdb::Connection;
use serde::Serialize;
use serde_json::to_writer_pretty;
use std::collections::HashMap;

use crate::sketch::sketching::*;

/// Arguments for the extract_sketches subcommand
#[derive(Args, Debug)]
pub struct ExtractSketchesArgs {
    /// Path to DuckDB database
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,

    /// Output JSON file path
    #[arg(long, required = true, value_name = "OUTPUT_JSON")]
    pub output: String,
}

/// Small container for export
#[derive(Serialize)]
struct Export {
    metadata: HashMap<String, String>,
    sketches: Vec<Sketch>,
}

/// Load metadata table into a HashMap<String,String>
fn load_metadata_map(conn: &Connection) -> Result<HashMap<String, String>> {
    // Query column names for metadata table
    let mut cols_stmt = conn.prepare(
        "SELECT column_name FROM information_schema.columns \
         WHERE table_name = 'metadata' ORDER BY ordinal_position"
    )?;

    let mut cols_rows = cols_stmt.query([])?;
    let mut cols: Vec<String> = Vec::new();
    while let Some(row) = cols_rows.next()? {
        let col_name: String = row.get(0)?;
        cols.push(col_name);
    }

    if cols.is_empty() {
        return Ok(HashMap::new());
    }

    let select_sql = format!("SELECT {} FROM metadata LIMIT 1", cols.join(", "));
    let mut row_stmt = conn.prepare(&select_sql)?;
    let map = row_stmt.query_row([], |r| {
        let mut m = HashMap::new();
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
            m.insert(col.clone(), val_str);
        }
        Ok(m)
    })?;

    Ok(map)
}

/// Export all sketches from the `sketches` table into a JSON file.
pub fn extract_sketches_command(args: &ExtractSketchesArgs) -> Result<()> {
    println!("=== Starting extract_sketches command ===");
    println!("Database: {}", args.db);
    println!("Output JSON: {}", args.output);

    if !Path::new(&args.db).exists() {
        return Err(anyhow::anyhow!("Database not found: {}", args.db));
    }

    // Open connection
    let conn = Connection::open(&args.db)
        .with_context(|| format!("Error opening database: {}", args.db))?;

    // Load metadata
    let metadata_map = load_metadata_map(&conn)
        .with_context(|| "Error loading metadata table")?;

    // Load SketchManager from DB
    let sketch_manager = load_sketch_manager_from_db(&conn, 21, 1000)
        .with_context(|| "Error loading SketchManager from DB")?;

    println!("✓ SketchManager loaded. Sketch count: {}", sketch_manager.count());

    // Extract sketches
    let sketches = sketch_manager.iter_sketches();

    // Prepare JSON output
    let file = File::create(&args.output)
        .with_context(|| format!("Error creating output file: {}", args.output))?;
    let writer = BufWriter::new(file);

    let export = Export {
        metadata: metadata_map,
        sketches,
    };

    // Serialize as pretty JSON. Sketch struct derives Serialize so `sketch` field (Vec<u64>) will be exported as array.
    to_writer_pretty(writer, &export)
        .with_context(|| format!("Error writing JSON to {}", args.output))?;

    println!("✓ Exported sketches and metadata to {}", args.output);
    Ok(())
}
