use clap::{Parser, Subcommand};
use std::process;
use anyhow::{Context, Result};

mod db;
mod sketch;
mod graph;
mod commands;

use commands::{init_command, update_command, verify_init_files};
use commands::update::UpdateArgs;

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
        #[arg(value_name = "TOML_FILE")]
        toml_file: String,

        /// Verificar la integridad de los archivos después de la creación
        #[arg(long)]
        verify: bool,
    },
    /// Actualiza la base de datos existente
    Update(UpdateArgs),  // <- Cambia aquí: usa UpdateArgs directamente
    /// Clasifica muestras usando la base de datos (por implementar)
    Classify {
        #[arg(value_name = "ACRONYM")]
        acronym: String,

        #[arg(value_name = "FASTA_FILE")]
        fasta_file: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Init { toml_file, verify } => {
            init_command(toml_file)
                .with_context(|| format!("Error ejecutando init para '{}'", toml_file))?;

            if *verify {
                let acronym = extract_acronym_from_toml(toml_file)
                    .with_context(|| format!("Error extrayendo acrónimo de {}", toml_file))?;
                verify_init_files(&acronym)
                    .with_context(|| format!("Error verificando init para '{}'", acronym))?;
                println!("✓ Verificación completada exitosamente");
            }
        }
        Commands::Update(args) => {
            // Llama a tu función update_command pasando args
            update_command(args)
                .with_context(|| "Error ejecutando update")?;
        }
        Commands::Classify { acronym, fasta_file } => {
            println!("Comando classify por implementar");
            println!("  - Proyecto: {}", acronym);
            println!("  - Archivo FASTA: {}", fasta_file);
            println!("Archivos esperados:");
            println!("  - {}.db", acronym);
            println!("  - {}_sketches.bin", acronym);
        }
    }

    Ok(())
}

/// Función auxiliar para extraer el acronym de un archivo TOML
fn extract_acronym_from_toml(toml_path: &str) -> Result<String> {
    use std::fs;
    use crate::db::db::MetadataConfig;

    let toml_content = fs::read_to_string(toml_path)
        .with_context(|| format!("No se pudo leer {}", toml_path))?;
    let config: MetadataConfig = toml::from_str(&toml_content)
        .with_context(|| format!("No se pudo parsear TOML {}", toml_path))?;
    Ok(config.acronym)
}
