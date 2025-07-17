use std::path::Path;

// Importar tu módulo
use bactaxid::db::db::*;

#[test]
fn test_full_database_workflow() {
    // Test con archivo temporal
    let db_path = "test_db.duckdb";
    
    {
        let db = DuckDb::new(db_path).expect("No se pudo crear la base de datos");
        
        // Inicializar datos de prueba
        let distances = vec![500, 1000, 1500, 2000];
        
        // Crear todas las tablas
        db.init_edges_table().expect("Error creando tabla edges");
        db.init_levels_table(distances).expect("Error creando tabla levels");
        db.init_code_table().expect("Error creando tabla code");
        db.init_code_full_table().expect("Error creando tabla code_full");
        db.init_code_state_table().expect("Error creando tabla code_state");
        
        // Verificar persistencia
        let conn = db.connection();
        let mut stmt = conn.prepare("SELECT COUNT(*) FROM levels").expect("Error preparando consulta");
        let mut rows = stmt.query([]).expect("Error ejecutando consulta");
        
        if let Some(row) = rows.next().expect("Error leyendo fila") {
            let count: i64 = row.get(0).expect("Error obteniendo count");
            assert_eq!(count, 4);
        }
    }
    
    // Limpiar archivo de prueba
    if Path::new(db_path).exists() {
        std::fs::remove_file(db_path).expect("Error eliminando archivo de prueba");
    }
}
