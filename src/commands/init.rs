

use crate::db::db::{DuckDb, MetadataConfig};

use anyhow::{Context, Result};
use duckdb::Result as DuckResult;
use std::fs;
use std::path::Path;



/// Verifica que el TOML existe, lo lee, parsea, crea e inicializa la base de datos.
/// Devuelve la instancia de DuckDb lista para usar.
pub fn verify_init_files(toml_path: &str) -> Result<DuckDb> {
    // 1. Verificar que el archivo TOML existe
    Path::new(toml_path)
        .exists()
        .then(|| ())
        .ok_or_else(|| anyhow::anyhow!("El archivo TOML no existe: {}", toml_path))?;

    // 2. Leer y parsear el archivo TOML
    let toml_content = fs::read_to_string(toml_path)
        .with_context(|| format!("Error leyendo TOML {}", toml_path))?;
    let config: MetadataConfig = toml::from_str(&toml_content)
        .with_context(|| format!("Error parseando TOML {}", toml_path))?;

    // 3. Crear la base de datos con el acronym
    let db_name = format!("{}.db", config.acronym);
    let db = DuckDb::new(&db_name)
        .with_context(|| format!("Error creando DuckDb {}", db_name))?;

    // 4. Inicializar todas las tablas según el config
    db.init_all_tables_from_config(&config)
        .with_context(|| format!("Error inicializando tablas en {}", db_name))?;

    // 5. Verificar que el archivo de la base de datos fue creado
    Path::new(&db_name)
        .exists()
        .then(|| ())
        .ok_or_else(|| anyhow::anyhow!("Error: no se creó la base de datos {}", db_name))?;

    Ok(db)
}

/// Función principal del comando init
pub fn init_command(toml_path: &str) -> Result<()> {
    println!("=== Iniciando comando init ===");

    // Usa la función de verificación e inicialización
    let db = verify_init_files(toml_path)?;

    // Mostrar información al usuario
    // Justo después de llamar a verify_init_files:
    let db_name = format!("{}.db", toml::from_str::<MetadataConfig>(&fs::read_to_string(toml_path)?).unwrap().acronym);

    let db_size = fs::metadata(&db_name)
        .with_context(|| format!("Error obteniendo tamaño de {}", db_name))?
        .len();

    println!("✓ Base de datos creada y cargada: {}", db_name);
    println!("  - Tamaño en disco: {} bytes", db_size);
    println!("=== Comando init completado exitosamente ===");

    Ok(())
}
