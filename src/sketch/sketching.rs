use needletail::parser::parse_fastx_file;
use needletail::sequence::Sequence;
use nthash::NtHashIterator;
use hashbrown::HashMap;
use bincode;
use serde::*;
use std::fs::File;
use std::io::{BufWriter,BufReader};
use rayon::prelude::*;
use anyhow::{Result, Context, anyhow, bail};
use duckdb::Connection;

use std::{
    error::Error,
    path::Path,
};


/// Estructura que representa un sketch de una secuencia biológica
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sketch {
    /// El sketch generado (vector de hashes mínimos)
    pub sketch: Vec<u64>,
    /// Tamaño de los k-mers utilizados
    pub kmer_size: usize,
    /// Tamaño del sketch (número de bins)
    pub sketch_size: usize,
    /// Nombre del archivo sin extensión
    pub name: String,
}

impl Sketch {
    /// Crea un nuevo sketch a partir de un archivo FASTA
    pub fn new<P: AsRef<Path>>(
        fasta_path: P,
        kmer_size: usize,
        sketch_size: usize,
    ) -> Result<Self> {
        // Generar el sketch usando la función binwise_minhash
        let sketch = binwise_minhash(&fasta_path, kmer_size, sketch_size)?;
        
        // Extraer el nombre del archivo sin extensión
        let name = fasta_path
            .as_ref()
            .file_stem()
            .and_then(|stem| stem.to_str())
            .unwrap_or("unknown")
            .to_string();

        Ok(Sketch {
            sketch,
            kmer_size,
            sketch_size,
            name,
        })
    }

    /// Calcula la distancia Mash entre este sketch y otro
    pub fn distance(&self, other: &Sketch) -> f64 {
        jaccard_index(&self.sketch, &other.sketch, self.kmer_size)
    }

    /// Verifica si dos sketches son compatibles para comparación
    pub fn is_compatible(&self, other: &Sketch) -> bool {
        self.kmer_size == other.kmer_size && self.sketch_size == other.sketch_size
    }

    /// Obtiene información básica del sketch
    pub fn info(&self) -> String {
        format!(
            "Sketch: {} (k={}, size={}, bins_filled={})",
            self.name,
            self.kmer_size,
            self.sketch_size,
            self.sketch.iter().filter(|&&x| x != u64::MAX).count()
        )
    }
}

/// Contenedor para gestionar múltiples sketches
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SketchManager {
    sketches: HashMap<String, Sketch>,
    /// Parámetros por defecto para nuevos sketches
    pub default_kmer_size: usize,
    pub default_sketch_size: usize,
}

impl SketchManager {
    /// 1. Inicializa un HashMap vacío para almacenar sketches
    pub fn new(default_kmer_size: usize, default_sketch_size: usize) -> Self {
        SketchManager {
            sketches: HashMap::new(),
            default_kmer_size,
            default_sketch_size,
        }
    }

    pub fn length(&self) -> usize {
        self.sketches.len()
    }
    /// 2. Añade un nuevo sketch al HashMap
    pub fn add_sketch(&mut self, sketch: Sketch) -> Result<()> {
        let name = sketch.name.clone();
        
        // ✅ SOLUCIÓN: Usar bail! en lugar de format!().into()
        if self.sketches.contains_key(&name) {
            bail!("Sketch con nombre '{}' ya existe", name);
        }

        self.sketches.insert(name, sketch);
        Ok(())
    }

    /// Añade un sketch desde un archivo FASTA
    pub fn add_sketch_from_file<P: AsRef<Path>>(
        &mut self,
        fasta_path: P,
        
    ) -> Result<()> {
            let k = self.default_kmer_size;
            let s = self.default_sketch_size;
        
        let sketch = Sketch::new(fasta_path, k, s)?;
        self.add_sketch(sketch)
    }

    /// Obtiene un sketch por nombre
    pub fn get_sketch(&self, name: &str) -> Option<&Sketch> {
        self.sketches.get(name)
    }

    /// Elimina un sketch por nombre
    pub fn remove_sketch(&mut self, name: &str) -> Option<Sketch> {
        self.sketches.remove(name)
    }

    /// Lista todos los nombres de sketches disponibles
    pub fn list_sketches(&self) -> Vec<String> {
        self.sketches.keys().cloned().collect()
    }

    /// Obtiene el número de sketches almacenados
    pub fn count(&self) -> usize {
        self.sketches.len()
    }

    pub fn contains(&self, name: &str) -> bool {
        self.sketches.contains_key(name)
    }

    /// Calcula la distancia entre dos sketches por nombre
    pub fn distance_between(&self, name1: &str, name2: &str) -> Result<f64> {
        let sketch1 = self.sketches.get(name1)
            .ok_or_else(|| anyhow!("Sketch '{}' no encontrado", name1))?;  // ← así sí

        let sketch2 = self.sketches.get(name2)
            .ok_or_else(|| anyhow!("Sketch '{}' no encontrado", name2))?;

        if !sketch1.is_compatible(sketch2) {
            anyhow::bail!("Los sketches no son compatibles (diferentes k-mer o sketch size)");
        }

        Ok(sketch1.distance(sketch2))
    }

    /// 3. Guarda el HashMap en disco usando bincode
    // pub fn save_to_disk<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
    //     // Crear el archivo de salida
    //     let file = File::create(file_path)?;
    //     let writer = BufWriter::new(file);
        
    //     // Serializar usando bincode
    //     bincode::serialize_into(writer, self)?;
        
    //     Ok(())
    // }

    // /// 4. Carga el HashMap desde disco usando bincode
    // pub fn load_from_disk<P: AsRef<Path>>(file_path: P) -> Result<Self> {
    //     // Abrir el archivo de entrada
    //     let file = File::open(file_path)?;
    //     let reader = BufReader::new(file);
        
    //     // Deserializar usando bincode
    //     let sketch_manager: SketchManager = bincode::deserialize_from(reader)?;
        
    //     Ok(sketch_manager)
    // }

    /// Función auxiliar para verificar si existe un archivo
    pub fn file_exists<P: AsRef<Path>>(file_path: P) -> bool {
        file_path.as_ref().exists()
    }


}


pub fn load_sketch_manager_from_db(
    conn: &Connection,
    default_kmer_size: usize,
    default_sketch_size: usize
) -> Result<SketchManager> {
    // 1. Obtener todos los datos serializados de forma secuencial
    let mut stmt = conn.prepare("SELECT sample, sketch FROM sketches")?;
    let rows = stmt.query_map([], |row| {
        let sample_id: String = row.get(0)?;
        let sketch_data: Vec<u8> = row.get(1)?;
        Ok((sample_id, sketch_data))
    })?;

    // 2. Recoger en un Vec para paralelizar
    let raw_data: Result<Vec<_>, _> = rows.collect();
    let raw_data = raw_data?;

    // Si no hay filas, crear un SketchManager vacío
    if raw_data.is_empty() {
        return Ok(SketchManager {
            sketches: HashMap::new(),
            default_kmer_size,
            default_sketch_size,
        });
    }

    // 3. Paralelizar la deserialización y construcción del HashMap
    let sketches: HashMap<String, Sketch> = raw_data
        .into_par_iter()
        .map(|(sample_id, sketch_data)| {
            let sketch: Sketch = bincode::deserialize(&sketch_data)
                .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;
            Ok((sample_id, sketch))
        })
        .collect::<Result<HashMap<_, _>, duckdb::Error>>()?;

    // 4. Construir el SketchManager
    Ok(SketchManager {
        sketches,
        default_kmer_size,
        default_sketch_size,
    })
}

/// Construye un b
/// - Usa binning directo en lugar de heap O(log S).
/// - Filtra duplicados por secuencia.
/// - Itera en streaming sin colectar todos los hashes.
/// - Derivado de BinDash: https://github.com/zhaoxiaofei/bindash
pub fn binwise_minhash<P: AsRef<Path>>(
    fasta_path: P,
    k: usize,
    sketch_size: usize,
) -> Result<Vec<u64>> {
    // Número de bins = sketch_size; tamaño del espacio de hash (u64::MAX+1)
    let bins = sketch_size;
    let bin_size = (u64::MAX / bins as u64).saturating_add(1);
    // Array de mínimos por bin
    let mut signs = vec![u64::MAX; bins];

    let mut reader = parse_fastx_file(fasta_path)?;
    while let Some(record) = reader.next() {
        let seqrec = record?;
        let seq_bytes = seqrec.normalize(false);

        // Rolling-hash streaming
        let mut iter = NtHashIterator::new(&seq_bytes, k)?;
        while let Some(hash) = iter.next() {
            // Mapear hash al bin correspondiente
            let idx = (hash / bin_size) as usize;
            // Mantener el mínimo en O(1)
            if hash < signs[idx] {
                signs[idx] = hash;
            }
        }
    }
    Ok(signs)
}


pub fn jaccard_index(sketch1: &[u64], sketch2: &[u64], kmer_size: usize) -> f64 {
    assert_eq!(sketch1.len(), sketch2.len(), "Sketches deben tener el mismo tamaño.");

    let matches = sketch1.iter()
        .zip(sketch2.iter())
        .filter(|(a, b)| a == b)
        .count();

    let jaccard_value = matches as f64 / sketch1.len() as f64;

    if jaccard_value == 0.0 {
        return 1.0; // Distancia máxima
    }

      // Calcular distancia Mash: D = -1/k * ln(2j/(1+j))
    let mash_distance = -((2.0 * jaccard_value) / (1.0 + jaccard_value)).ln() / kmer_size as f64;
    
    // Convertir distancia Mash a ANI
    let ani = 1.0 - mash_distance;
    
    // Asegurar que ANI esté en el rango [0, 1]
    ani.max(0.0).min(1.0)
}
/// Compara los sketches de query_manager contra un subconjunto específico de sketches 
/// en reference_manager, devolviendo solo aquellas comparaciones con distancia < max_dist.
///
/// # Argumentos
/// - `query_manager`: SketchManager con los sketches fuente (query).
/// - `reference_manager`: SketchManager con los sketches de referencia.
/// - `subset_ids`: Lista de IDs específicos para comparar del reference_manager.
/// - `max_dist`: Distancia máxima. Solo se retornan distancias < max_dist.
///
/// # Retorna
/// Vector de tuplas (nombre_query, nombre_referencia, distancia) donde distancia < max_dist.
pub fn pw_one_to_many(
    query_manager: &SketchManager,
    reference_manager: &SketchManager,
    subset_ids: &[String],
    min_dist: f64,
) -> Vec<(String, String, f64)> {
    query_manager.sketches.par_iter().flat_map(|(source_name, source_sketch)| {
        subset_ids.par_iter().filter_map(move |target_id| {
            reference_manager.get_sketch(target_id).and_then(|target_sketch| {
                let dist = source_sketch.distance(target_sketch);
                if dist > min_dist {
                    Some((source_name.clone(), target_id.clone(), dist))
                } else {
                    None
                }
            })
        }).collect::<Vec<_>>()
    }).collect()
}

/// Compara todos los sketches en `query_manager` contra todos los sketches en `reference_manager`
/// de forma paralela, devolviendo solo aquellas comparaciones con distancia < max_dist.
///
/// # Argumentos
/// - `query_manager`: Referencia al SketchManager que contiene los sketches fuente (query).
/// - `reference_manager`: Referencia al SketchManager que contiene los sketches de referencia (target).
/// - `max_dist`: Distancia máxima. Solo se retornan distancias < max_dist.
///
/// # Retorna
/// Un vector de tuplas (nombre_source, nombre_target, distancia) donde distancia < max_dist.
pub fn pw_one_to_all(
    query_manager: &SketchManager,
    reference_manager: &SketchManager,
    max_dist: f64,
) -> Vec<(String, String, f64)> {
    query_manager.sketches.par_iter().flat_map(|(source_name, source_sketch)| {
        reference_manager.sketches.par_iter().filter_map(move |(target_name, target_sketch)| {
            let dist = source_sketch.distance(target_sketch);
            if dist < max_dist {
                Some((source_name.clone(), target_name.clone(), dist))
            } else {
                None
            }
        })
    }).collect()
}
