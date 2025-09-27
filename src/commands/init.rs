

use crate::db::db::{DuckDb, MetadataConfig};

use anyhow::{Context, Result};
use duckdb::Result as DuckResult;
use std::fs;
use std::path::Path;


/// Checks that the TOML exists, reads it, parses it, creates and initializes the database.
/// Returns the DuckDb instance ready to use.
pub fn verify_init_files(toml_path: &str) -> Result<DuckDb> {
    // 1. Check that the TOML file exists
    Path::new(toml_path)
        .exists()
        .then(|| ())
        .ok_or_else(|| anyhow::anyhow!("TOML file does not exist: {}", toml_path))?;

    // 2. Read and parse the TOML file
    let toml_content = fs::read_to_string(toml_path)
        .with_context(|| format!("Error reading TOML {}", toml_path))?;
    let config: MetadataConfig = toml::from_str(&toml_content)
        .with_context(|| format!("Error parsing TOML {}", toml_path))?;

    // 3. Create the database with the acronym
    let db_name = format!("{}.db", config.acronym);
    let db = DuckDb::new(&db_name)
        .with_context(|| format!("Error creating DuckDb {}", db_name))?;

    // 4. Initialize all tables according to config
    db.init_all_tables_from_config(&config)
        .with_context(|| format!("Error initializing tables in {}", db_name))?;

    // 5. Check that the database file was created
    Path::new(&db_name)
        .exists()
        .then(|| ())
        .ok_or_else(|| anyhow::anyhow!("Error: database file {} was not created", db_name))?;

    Ok(db)
}

/// Main function for the init command
pub fn init_command(toml_path: &str) -> Result<()> {
    println!("=== Starting init command ===");

    // Use the verification and initialization function
    let db = verify_init_files(toml_path)?;

    // Show information to the user
    // Right after calling verify_init_files:
    let db_name = format!("{}.db", toml::from_str::<MetadataConfig>(&fs::read_to_string(toml_path)?).unwrap().acronym);

    let db_size = fs::metadata(&db_name)
        .with_context(|| format!("Error getting size of {}", db_name))?
        .len();

    println!("✓ Database created and loaded: {}", db_name);
    println!("  - Disk size: {} bytes", db_size);
    println!("=== Init command successfully completed ===");

    Ok(())
}
