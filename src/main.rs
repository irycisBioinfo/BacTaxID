use clap::{Parser, Subcommand};
use anyhow::{Context, Result};

mod db;
mod sketch;
mod graph;
mod commands;

use commands::{init_command, update_command, verify_init_files};  // ADDED
use commands::update::UpdateArgs;
use commands::classify::ClassifyArgs;
use commands::distance::DistanceArgs;  // ADDED
use commands::extract_sketches::ExtractSketchesArgs;

#[derive(Parser)]
#[command(name = "bactaxid")]
#[command(about = "Bacterial taxonomic analysis tool")]
#[command(version = "1.0")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Initializes a new database and SketchManager from a TOML file
    Init {
        #[arg(value_name = "TOML_FILE")]
        toml_file: String,
        
        /// Verify the integrity of files after creation
        #[arg(long)]
        verify: bool,
    },
    
    /// Updates the existing database
    Update(UpdateArgs),
    
    /// Classifies sequences against the reference database (read-only)
    Classify(ClassifyArgs),
    
    /// Calculates pairwise distances between sequences (read-only)
    Distance(DistanceArgs),  // ADDED
    /// Export all sketches from DB into JSON
    ExtractSketches(ExtractSketchesArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    
    match &cli.command {
        Commands::Init { toml_file, verify } => {
            init_command(toml_file)
                .with_context(|| format!("Error running init for '{}'", toml_file))?;
            
            if *verify {
                let acronym = extract_acronym_from_toml(toml_file)
                    .with_context(|| format!("Error extracting acronym from {}", toml_file))?;
                verify_init_files(&acronym)
                    .with_context(|| format!("Error verifying init for '{}'", acronym))?;
                println!("✓ Verification successfully completed");
            }
        }
        
        Commands::Update(args) => {
            update_command(args)
                .with_context(|| "Error running update")?;
        }
        
        Commands::Classify(args) => {
            use commands::classify::classify_command;
            classify_command(args)
                .with_context(|| "Error running classify")?;
        }
        
        Commands::Distance(args) => {  // ADDED
            use commands::distance::distance_command;
            distance_command(args)
                .with_context(|| "Error running distance")?;
        }
        Commands::ExtractSketches(args) => {
            use commands::extract_sketches::extract_sketches_command;
            extract_sketches_command(args)
                .with_context(|| "Error running extract_sketches")?;
        }
    }
    
    Ok(())
}

/// Helper function to extract acronym from a TOML file
fn extract_acronym_from_toml(toml_path: &str) -> Result<String> {
    use std::fs;
    use crate::db::db::MetadataConfig;
    
    let toml_content = fs::read_to_string(toml_path)
        .with_context(|| format!("Could not read {}", toml_path))?;
    
    let config: MetadataConfig = toml::from_str(&toml_content)
        .with_context(|| format!("Could not parse TOML {}", toml_path))?;
    
    Ok(config.acronym)
}
