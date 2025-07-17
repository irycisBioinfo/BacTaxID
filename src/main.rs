use clap::{Parser, Subcommand};
use std::process;

mod db;
mod sketch;
mod graph;
mod commands;

use commands::{init_command, verify_init_files};

#[derive(Parser)]
#[command(name = "bactaxid")]
#[command(about = "Herramienta de análisis taxonómico bacteriano")]
#[command(version = "1.0")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Inicializa una nueva base de datos y SketchManager desde un archivo TOML
    Init {
        /// Ruta al archivo TOML de configuración
        #[arg(value_name = "TOML_FILE")]
        toml_file: String,
        
        /// Verificar la integridad de los archivos después de la creación
        #[arg(long)]
        verify: bool,
    },
    /// Actualiza la base de datos existente (por implementar)
    Update {
        /// Acronym del proyecto a actualizar
        #[arg(value_name = "ACRONYM")]
        acronym: String,
    },
    /// Clasifica muestras usando la base de datos (por implementar)
    Classify {
        /// Acronym del proyecto a usar para clasificación
        #[arg(value_name = "ACRONYM")]
        acronym: String,
        
        /// Archivo FASTA con las muestras a clasificar
        #[arg(value_name = "FASTA_FILE")]
        fasta_file: String,
    },
}

fn main() {
    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::Init { toml_file, verify } => {
            match init_command(toml_file) {
                Ok(_) => {
                    if *verify {
                        // Extraer el acronym del archivo TOML para verificación
                        match extract_acronym_from_toml(toml_file) {
                            Ok(acronym) => {
                                match verify_init_files(&acronym) {
                                    Ok(_) => {
                                        println!("✓ Verificación completada exitosamente");
                                        Ok(())
                                    }
                                    Err(e) => Err(e),
                                }
                            }
                            Err(e) => Err(e),
                        }
                    } else {
                        Ok(())
                    }
                }
                Err(e) => Err(e),
            }
        }
        Commands::Update { acronym } => {
            println!("Comando update por implementar para: {}", acronym);
            println!("Archivos esperados:");
            println!("  - {}.db", acronym);
            println!("  - {}_sketches.bin", acronym);
            Ok(())
        }
        Commands::Classify { acronym, fasta_file } => {
            println!("Comando classify por implementar");
            println!("  - Proyecto: {}", acronym);
            println!("  - Archivo FASTA: {}", fasta_file);
            println!("Archivos esperados:");
            println!("  - {}.db", acronym);
            println!("  - {}_sketches.bin", acronym);
            Ok(())
        }
    };

    if let Err(e) = result {
        eprintln!("Error: {}", e);
        process::exit(1);
    }
}

/// Función auxiliar para extraer el acronym de un archivo TOML
fn extract_acronym_from_toml(toml_path: &str) -> Result<String, Box<dyn std::error::Error>> {
    use std::fs;
    use crate::db::db::MetadataConfig;
    
    let toml_content = fs::read_to_string(toml_path)?;
    let config: MetadataConfig = toml::from_str(&toml_content)?;
    Ok(config.acronym)
}
