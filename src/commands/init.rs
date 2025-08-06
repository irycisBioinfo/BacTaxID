use duckdb::Result;
use std::path::Path;
use std::fs;

use crate::db::db::{DuckDb, MetadataConfig};
use crate::sketch::sketching::SketchManager;

/// Función principal del comando init
pub fn init_command(toml_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Iniciando comando init ===");
    
    // 1. Validar que el archivo TOML existe
    if !Path::new(toml_path).exists() {
        return Err(format!("El archivo TOML no existe: {}", toml_path).into());
    }
    
    // 2. Leer el archivo TOML para obtener los parámetros
    let toml_content = fs::read_to_string(toml_path)?;
    let config: MetadataConfig = toml::from_str(&toml_content)?;
    
    println!("✓ Archivo TOML leído correctamente");
    println!("  - Genus: {}", config.genus);
    println!("  - Acronym: {}", config.acronym);
    println!("  - K-mer size: {}", config.kmer_size);
    println!("  - Sketch size: {}", config.sketch_size);
    println!("  - Click size: {}", config.click_size);
    println!("  - Click threshold: {}", config.click_threshold);
    println!("  - Reference Size:{}", config.reference_size);
    
    // 3. Crear la base de datos usando el acronym (no genus)
    let db_name = format!("{}.db", config.acronym);
    let db = DuckDb::new(&db_name)?;
    
    // 4. Inicializar todas las tablas usando los datos del TOML
    db.init_all_tables_from_config(&config)
        .map_err(|e| format!("Error inicializando tablas: {}", e))?;
    
    println!("✓ Base de datos creada: {}", db_name);
    println!("✓ Todas las tablas inicializadas correctamente");
    
    // 5. Crear un SketchManager vacío usando los parámetros del TOML
    let sketch_manager = SketchManager::new(
        config.kmer_size as usize,
        config.sketch_size as usize
    );
    
    println!("✓ SketchManager vacío creado");
    println!("  - Default k-mer size: {}", config.kmer_size);
    println!("  - Default sketch size: {}", config.sketch_size);
    
    // 6. Guardar el SketchManager en disco usando el acronym
    let sketch_file = format!("{}_sketches.bin", config.acronym);
    sketch_manager.save_to_disk(&sketch_file)?;
    
    println!("✓ SketchManager guardado en: {}", sketch_file);
    
    // 7. Verificar que los archivos se crearon correctamente
    if Path::new(&db_name).exists() && Path::new(&sketch_file).exists() {
        println!("=== Comando init completado exitosamente ===");
        println!("Archivos creados:");
        println!("  - Base de datos: {} ({} bytes)", db_name, fs::metadata(&db_name)?.len());
        println!("  - SketchManager: {} ({} bytes)", sketch_file, fs::metadata(&sketch_file)?.len());
    } else {
        return Err("Error: Los archivos no se crearon correctamente".into());
    }
    
    Ok(())
}

/// Función auxiliar para verificar la integridad de los archivos creados
pub fn verify_init_files(acronym: &str) -> Result<(), Box<dyn std::error::Error>> {
    let db_name = format!("{}.db", acronym);
    let sketch_file = format!("{}_sketches.bin", acronym);
    
    // Verificar que la base de datos existe y es válida
    if !Path::new(&db_name).exists() {
        return Err(format!("Base de datos no encontrada: {}", db_name).into());
    }
    
    // Intentar conectar a la base de datos
    let db = DuckDb::new(&db_name)?;
    
    // Verificar que las tablas existen
    let tables_query = "SELECT table_name FROM information_schema.tables WHERE table_name IN ('edges', 'levels', 'code', 'code_full', 'code_state', 'metadata') ORDER BY table_name";
    let mut stmt = db.connection().prepare(tables_query)?;
    let mut rows = stmt.query([])?;
    
    let mut tables = Vec::new();
    while let Some(row) = rows.next()? {
        let table_name: String = row.get(0)?;
        tables.push(table_name);
    }
    
    let expected_tables = vec!["code", "code_full", "code_state", "edges", "levels", "metadata"];
    if tables != expected_tables {
        return Err(format!("Tablas faltantes en la base de datos. Encontradas: {:?}, Esperadas: {:?}", tables, expected_tables).into());
    }
    
    // Verificar que el SketchManager existe y es válido
    if !Path::new(&sketch_file).exists() {
        return Err(format!("SketchManager no encontrado: {}", sketch_file).into());
    }
    
    // Intentar cargar el SketchManager
    let _sketch_manager = SketchManager::load_from_disk(&sketch_file)?;
    
    println!("✓ Verificación completada - Todos los archivos son válidos");
    Ok(())
}
