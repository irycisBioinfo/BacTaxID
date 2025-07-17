// src/main.rs
use bactaxid::db::db::*;  // Reemplaza 'tu_nombre_proyecto' con el nombre real

fn main() -> duckdb::Result<()> {
    let db = DuckDb::new("mydb.duckdb")?;
    
    // Inicializar tablas
    let distances = vec![1000, 2500, 3200, 4800];
    db.init_levels_table(distances)?;
    db.init_edges_table()?;
    db.init_code_table()?;
    db.init_code_full_table()?;
    db.init_code_state_table()?;
    
    println!("Base de datos inicializada correctamente.");
    Ok(())
}