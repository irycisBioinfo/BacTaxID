use needletail::parser::parse_fastx_file;
use needletail::sequence::Sequence;
use nthash::NtHashIterator;
use hashbrown::HashMap;
use serde::*;
use rayon::prelude::*;
use anyhow::{Result, anyhow, bail};
use duckdb::Connection;
use twox_hash::XxHash64; 
use std::hash::Hasher; 



use std::path::Path;


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
    /// xxHash64 signature of the sketch vector
    pub signature: u64,
}

// Alias for compatibility with db.rs
// This allows db.rs to use sketch.hashes while keeping internal field as sketch.sketch
impl Sketch {
    pub fn hashes(&self) -> &[u64] {
        &self.sketch
    }
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


        // Calculate signature from sketch vector
        let signature = compute_sketch_signature(&sketch);


        Ok(Sketch {
            sketch,
            kmer_size,
            sketch_size,
            name,
            signature,
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
            "Sketch: {} (k={}, size={}, bins_filled={}, signature={})",
            self.name,
            self.kmer_size,
            self.sketch_size,
            self.sketch.iter().filter(|&&x| x != u64::MAX).count(),
            self.signature
        )
    }
}


/// Container to manage multiple sketches
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SketchManager {
    sketches: HashMap<u64, Sketch>,
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
        let signature = sketch.signature;


        // ✅ Check if signature already exists
        if self.sketches.contains_key(&signature) {
            bail!("Sketch with signature '{}' already exists", signature);
        }


        self.sketches.insert(signature, sketch);
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


    /// Gets a sketch by signature
    pub fn get_sketch(&self, signature: u64) -> Option<&Sketch> {
        self.sketches.get(&signature)
    }


    /// Gets a sketch by name (helper function)
    pub fn get_sketch_by_name(&self, name: &str) -> Option<&Sketch> {
        self.sketches.values().find(|s| s.name == name)
    }


    /// Removes a sketch by signature
    pub fn remove_sketch(&mut self, signature: u64) -> Option<Sketch> {
        self.sketches.remove(&signature)
    }


    /// Lists all available sketch names
    pub fn list_sketches(&self) -> Vec<String> {
        self.sketches.values().map(|s| s.name.clone()).collect()
    }


    /// Lists all available sketch signatures
    pub fn list_signatures(&self) -> Vec<u64> {
        self.sketches.keys().cloned().collect()
    }


    /// Gets the number of stored sketches
    pub fn count(&self) -> usize {
        self.sketches.len()
    }


    pub fn contains(&self, signature: u64) -> bool {
        self.sketches.contains_key(&signature)
    }


    pub fn contains_by_name(&self, name: &str) -> bool {
        self.sketches.values().any(|s| s.name == name)
    }


    /// Calculates the distance between two sketches by signature
    pub fn distance_between(&self, sig1: u64, sig2: u64) -> Result<f64> {
        let sketch1 = self.sketches.get(&sig1)
            .ok_or_else(|| anyhow!("Sketch with signature '{}' not found", sig1))?;


        let sketch2 = self.sketches.get(&sig2)
            .ok_or_else(|| anyhow!("Sketch with signature '{}' not found", sig2))?;


        if !sketch1.is_compatible(sketch2) {
            anyhow::bail!("Sketches are not compatible (different k-mer or sketch size)");
        }


        Ok(sketch1.distance(sketch2))
    }
}


/// Loads SketchManager from database using BLOB
pub fn load_sketch_manager_from_db(
    conn: &Connection,
    default_kmer_size: usize,
    default_sketch_size: usize
) -> Result<SketchManager> {
    // 1. Get all data sequentially - read BLOB directly
    let mut stmt = conn.prepare("SELECT sketch FROM sketches")?;
    let rows = stmt.query_map([], |row| {
        let sketch_blob: Vec<u8> = row.get(0)?;
        Ok(sketch_blob)
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

    // 3. Parallelize Sketch deserialization
    let sketches: HashMap<u64, Sketch> = raw_data
        .into_par_iter()
        .map(|sketch_blob| {
            // Deserialize Sketch directly with bincode
            let sketch: Sketch = bincode::deserialize(&sketch_blob)
                .map_err(|e| anyhow!("Bincode deserialization error: {}", e))?;
            
            Ok((sketch.signature, sketch))
        })
        .collect::<Result<HashMap<_, _>, anyhow::Error>>()?;

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


/// Computes xxHash64 signature of a sketch vector
/// This function creates a unique identifier for each sketch
fn compute_sketch_signature(sketch: &[u64]) -> u64 {
    let mut hasher = XxHash64::with_seed(0);
    
    // Hash each u64 in little-endian format for reproducibility
    for &value in sketch {
        hasher.write(&value.to_le_bytes());
    }
    
    hasher.finish()
}


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
    } else {
        // Compute Mash distance: D = -1/k * ln(2j/(1+j))
        let mash_distance = -((2.0 * jaccard_value) / (1.0 + jaccard_value)).ln() / kmer_size as f64;


        // Convert Mash distance to ANI
        let ani = 1.0 - mash_distance;


        // Make sure ANI is in [0, 1] range
        return ani.max(0.0).min(1.0);
    }
}


/// Compares the query_manager sketches against a specific subset of sketches 
/// in reference_manager, only returning those comparisons with distance > min_dist.
///
/// # Arguments
/// - `query_manager`: SketchManager with query/source sketches.
/// - `reference_manager`: SketchManager with reference sketches.
/// - `subset_sigs`: List of specific signatures to compare from reference_manager.
/// - `min_dist`: Minimum distance. Only distances > min_dist are returned.
///
/// # Returns
/// Vector of tuples (query_name, reference_name, distance) where distance > min_dist.
pub fn pw_one_to_many(
    query_manager: &SketchManager,
    reference_manager: &SketchManager,
    subset_sigs: &[u64],
    min_dist: f64,
) -> Vec<(u64, u64, f64)> {
    query_manager.sketches.par_iter().flat_map(|(query_sig, source_sketch)| {
        subset_sigs.par_iter().filter_map(move |target_sig| {
            reference_manager.get_sketch(*target_sig).and_then(|target_sketch| {
                let dist = source_sketch.distance(target_sketch);
                if dist > min_dist {
                    Some((*query_sig, *target_sig, dist))
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
    query_manager.sketches.par_iter().flat_map(|(_, source_sketch)| {
        reference_manager.sketches.par_iter().filter_map(move |(_, target_sketch)| {
            let dist = source_sketch.distance(target_sketch);
            if dist < max_dist {
                Some((source_sketch.name.clone(), target_sketch.name.clone(), dist))
            } else {
                None
            }
        })
    }).collect()
}
