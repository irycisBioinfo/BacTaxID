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


/// Structure that represents a sketch of a biological sequence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sketch {
    /// The generated sketch (vector of minimum hashes)
    pub sketch: Vec<u64>,
    /// Size of the k-mers used
    pub kmer_size: usize,
    /// Size of the sketch (number of bins)
    pub sketch_size: usize,
    /// Filename without extension
    pub name: String,
}

impl Sketch {
    /// Creates a new sketch from a FASTA file
    pub fn new<P: AsRef<Path>>(
        fasta_path: P,
        kmer_size: usize,
        sketch_size: usize,
    ) -> Result<Self> {
        // Generate the sketch using the binwise_minhash function
        let sketch = binwise_minhash(&fasta_path, kmer_size, sketch_size)?;

        // Extract the filename without extension
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

    /// Computes the Mash distance between this sketch and another
    pub fn distance(&self, other: &Sketch) -> f64 {
        jaccard_index(&self.sketch, &other.sketch, self.kmer_size)
    }

    /// Checks if two sketches are compatible for comparison
    pub fn is_compatible(&self, other: &Sketch) -> bool {
        self.kmer_size == other.kmer_size && self.sketch_size == other.sketch_size
    }

    /// Gets basic information about the sketch
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

/// Container to manage multiple sketches
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SketchManager {
    sketches: HashMap<String, Sketch>,
    /// Default parameters for new sketches
    pub default_kmer_size: usize,
    pub default_sketch_size: usize,
}

impl SketchManager {
    /// 1. Initializes an empty HashMap to store sketches
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
    /// 2. Adds a new sketch to the HashMap
    pub fn add_sketch(&mut self, sketch: Sketch) -> Result<()> {
        let name = sketch.name.clone();

        // ✅ SOLUTION: Use bail! instead of format!().into()
        if self.sketches.contains_key(&name) {
            bail!("Sketch with name '{}' already exists", name);
        }

        self.sketches.insert(name, sketch);
        Ok(())
    }

    /// Adds a sketch from a FASTA file
    pub fn add_sketch_from_file<P: AsRef<Path>>(
        &mut self,
        fasta_path: P,
    ) -> Result<()> {
        let k = self.default_kmer_size;
        let s = self.default_sketch_size;
        let sketch = Sketch::new(fasta_path, k, s)?;
        self.add_sketch(sketch)
    }

    /// Gets a sketch by name
    pub fn get_sketch(&self, name: &str) -> Option<&Sketch> {
        self.sketches.get(name)
    }

    /// Removes a sketch by name
    pub fn remove_sketch(&mut self, name: &str) -> Option<Sketch> {
        self.sketches.remove(name)
    }

    /// Lists all available sketch names
    pub fn list_sketches(&self) -> Vec<String> {
        self.sketches.keys().cloned().collect()
    }

    /// Gets the number of stored sketches
    pub fn count(&self) -> usize {
        self.sketches.len()
    }

    pub fn contains(&self, name: &str) -> bool {
        self.sketches.contains_key(name)
    }

    /// Calculates the distance between two sketches by name
    pub fn distance_between(&self, name1: &str, name2: &str) -> Result<f64> {
        let sketch1 = self.sketches.get(name1)
            .ok_or_else(|| anyhow!("Sketch '{}' not found", name1))?;

        let sketch2 = self.sketches.get(name2)
            .ok_or_else(|| anyhow!("Sketch '{}' not found", name2))?;

        if !sketch1.is_compatible(sketch2) {
            anyhow::bail!("Sketches are not compatible (different k-mer or sketch size)");
        }

        Ok(sketch1.distance(sketch2))
    }

    /// 3. Saves the HashMap to disk using bincode
    // pub fn save_to_disk<P: AsRef<Path>>(&self, file_path: P) -> Result<()> {
    //     // Create the output file
    //     let file = File::create(file_path)?;
    //     let writer = BufWriter::new(file);
    //     // Serialize using bincode
    //     bincode::serialize_into(writer, self)?;
    //     Ok(())
    // }

    // /// 4. Loads the HashMap from disk using bincode
    // pub fn load_from_disk<P: AsRef<Path>>(file_path: P) -> Result<Self> {
    //     // Open the input file
    //     let file = File::open(file_path)?;
    //     let reader = BufReader::new(file);
    //     // Deserialize using bincode
    //     let sketch_manager: SketchManager = bincode::deserialize_from(reader)?;
    //     Ok(sketch_manager)
    // }

    /// Helper function to check if a file exists
    pub fn file_exists<P: AsRef<Path>>(file_path: P) -> bool {
        file_path.as_ref().exists()
    }
}

pub fn load_sketch_manager_from_db(
    conn: &Connection,
    default_kmer_size: usize,
    default_sketch_size: usize
) -> Result<SketchManager> {
    // 1. Get all serialized data sequentially
    let mut stmt = conn.prepare("SELECT sample, sketch FROM sketches")?;
    let rows = stmt.query_map([], |row| {
        let sample_id: String = row.get(0)?;
        let sketch_data: Vec<u8> = row.get(1)?;
        Ok((sample_id, sketch_data))
    })?;

    // 2. Gather in a Vec for parallelization
    let raw_data: Result<Vec<_>, _> = rows.collect();
    let raw_data = raw_data?;

    // If no rows, create an empty SketchManager
    if raw_data.is_empty() {
        return Ok(SketchManager {
            sketches: HashMap::new(),
            default_kmer_size,
            default_sketch_size,
        });
    }

    // 3. Parallelize deserialization and HashMap construction
    let sketches: HashMap<String, Sketch> = raw_data
        .into_par_iter()
        .map(|(sample_id, sketch_data)| {
            let sketch: Sketch = bincode::deserialize(&sketch_data)
                .map_err(|e| duckdb::Error::ToSqlConversionFailure(Box::new(e)))?;
            Ok((sample_id, sketch))
        })
        .collect::<Result<HashMap<_, _>, duckdb::Error>>()?;

    // 4. Construct the SketchManager
    Ok(SketchManager {
        sketches,
        default_kmer_size,
        default_sketch_size,
    })
}

/// Constructs a binwise sketch
/// - Uses direct binning instead of heap O(log S).
/// - Filters duplicates by sequence.
/// - Iterates streaming without collecting all hashes.
/// - Derived from BinDash: [https://github.com/zhaoxiaofei/bindash](https://github.com/zhaoxiaofei/bindash)
pub fn binwise_minhash<P: AsRef<Path>>(
    fasta_path: P,
    k: usize,
    sketch_size: usize,
) -> Result<Vec<u64>> {
    // Number of bins = sketch_size; hash space size (u64::MAX+1)
    let bins = sketch_size;
    let bin_size = (u64::MAX / bins as u64).saturating_add(1);
    // Array of minimums per bin
    let mut signs = vec![u64::MAX; bins];

    let mut reader = parse_fastx_file(fasta_path)?;
    while let Some(record) = reader.next() {
        let seqrec = record?;
        let seq_bytes = seqrec.normalize(false);

        // Rolling-hash streaming
        let mut iter = NtHashIterator::new(&seq_bytes, k)?;
        while let Some(hash) = iter.next() {
            // Map hash to the corresponding bin
            let idx = (hash / bin_size) as usize;
            // Maintain the minimum in O(1)
            if hash < signs[idx] {
                signs[idx] = hash;
            }
        }
    }
    Ok(signs)
}

// This function actually returns ANI instead of Jaccard
/// Calculates the Jaccard index between two sketches and converts it to ANI.
/// - `sketch1` and `sketch2`: Vectors of minimum hashes
/// - `kmer_size`: Size of the k-mers used for the sketch.
/// - Returns the Jaccard index as a value between 0 and 1
/// - If the Jaccard index is 0, returns 0.0 (minimum similarity).
pub fn jaccard_index(sketch1: &[u64], sketch2: &[u64], kmer_size: usize) -> f64 {
    assert_eq!(sketch1.len(), sketch2.len(), "Sketches must have the same size.");

    let matches = sketch1.iter()
        .zip(sketch2.iter())
        .filter(|(a, b)| a == b)
        .count();

    let jaccard_value = matches as f64 / sketch1.len() as f64;

    if jaccard_value == 0.0 {
        return 0.0; // Minimum similarity
    } else{
        // Compute Mash distance: D = -1/k * ln(2j/(1+j))
        let mash_distance = -((2.0 * jaccard_value) / (1.0 + jaccard_value)).ln() / kmer_size as f64;

        // Convert Mash distance to ANI
        let ani = 1.0 - mash_distance;

        // Make sure ANI is in [0, 1] range
        return ani.max(0.0).min(1.0);
    }
}

/// Compares the query_manager sketches against a specific subset of sketches 
/// in reference_manager, only returning those comparisons with distance < max_dist.
///
/// # Arguments
/// - `query_manager`: SketchManager with query/source sketches.
/// - `reference_manager`: SketchManager with reference sketches.
/// - `subset_ids`: List of specific IDs to compare from reference_manager.
/// - `max_dist`: Maximum distance. Only distances < max_dist are returned.
///
/// # Returns
/// Vector of tuples (query_name, reference_name, distance) where distance < max_dist.
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

/// Compares all sketches in `query_manager` against all in `reference_manager`
/// in parallel, only returning comparisons with distance < max_dist.
///
/// # Arguments
/// - `query_manager`: Reference to the SketchManager containing source (query) sketches.
/// - `reference_manager`: Reference to the SketchManager containing reference (target) sketches.
/// - `max_dist`: Maximum distance. Only distances < max_dist are returned.
///
/// # Returns
/// Vector of tuples (source_name, target_name, distance) where distance < max_dist.
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
