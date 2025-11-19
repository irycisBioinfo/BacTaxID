use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use serde::{Deserialize, Serialize};
use anyhow::{Result, Context};

/// Represents a genomic sketch with MinHash values and metadata
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Sketch {
    pub name: String,
    pub signature: u64,        // Hash del sketch (identificador único)
    pub sketch: Vec<u64>,      // MinHash values
    pub kmer_size: usize,
    pub sketch_size: usize,
}

impl Sketch {
    /// Creates a new Sketch from a FASTA file
    pub fn new(path: &Path, kmer_size: usize, sketch_size: usize) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let name = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();

        let mut sequence = String::new();
        for line in reader.lines() {
            let line = line.context("Failed to read line from file")?;
            if !line.starts_with('>') {
                sequence.push_str(&line.trim());
            }
        }

        if sequence.is_empty() {
            anyhow::bail!("Empty sequence in file: {}", path.display());
        }

        // Generate k-mers and create MinHash sketch
        let sketch_values = create_minhash_sketch(&sequence, kmer_size, sketch_size)?;

        // Calculate signature as hash of the sketch
        let signature = calculate_sketch_signature(&sketch_values);

        Ok(Sketch {
            name,
            signature,
            sketch: sketch_values,
            kmer_size,
            sketch_size,
        })
    }

    /// Calculates distance (similarity) between two sketches using Jaccard index
    pub fn distance(&self, other: &Sketch) -> f64 {
        if self.kmer_size != other.kmer_size {
            eprintln!(
                "Warning: k-mer sizes differ ({} vs {})",
                self.kmer_size, other.kmer_size
            );
            return 0.0;
        }

        // Calculate Jaccard index
        let intersection = self
            .sketch
            .iter()
            .filter(|&hash| other.sketch.binary_search(hash).is_ok())
            .count();

        let union = self.sketch.len() + other.sketch.len() - intersection;

        if union == 0 {
            return 1.0;
        }

        let jaccard = intersection as f64 / union as f64;
        
        // Convert Jaccard to Mash distance and then to ANI
        let mash_distance = jaccard_to_mash(jaccard, self.kmer_size);
        
        // Return ANI (1 - mash_distance)
        1.0 - mash_distance
    }
}

/// Creates a MinHash sketch from a DNA sequence
fn create_minhash_sketch(sequence: &str, kmer_size: usize, sketch_size: usize) -> Result<Vec<u64>> {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    if sequence.len() < kmer_size {
        anyhow::bail!("Sequence length ({}) is shorter than k-mer size ({})", 
                     sequence.len(), kmer_size);
    }

    let sequence_upper = sequence.to_uppercase();
    let mut kmer_hashes = Vec::new();

    // Generate all k-mers and hash them
    for i in 0..=(sequence_upper.len() - kmer_size) {
        let kmer = &sequence_upper[i..i + kmer_size];
        
        // Skip k-mers with N's
        if kmer.contains('N') {
            continue;
        }

        // Hash the k-mer
        let mut hasher = DefaultHasher::new();
        kmer.hash(&mut hasher);
        let hash = hasher.finish();
        kmer_hashes.push(hash);
        
        // Also hash reverse complement
        let rev_comp = reverse_complement(kmer);
        let mut hasher_rc = DefaultHasher::new();
        rev_comp.hash(&mut hasher_rc);
        let hash_rc = hasher_rc.finish();
        kmer_hashes.push(hash_rc);
    }

    if kmer_hashes.is_empty() {
        anyhow::bail!("No valid k-mers generated from sequence");
    }

    // Sort and take the smallest sketch_size hashes (MinHash)
    kmer_hashes.sort_unstable();
    kmer_hashes.dedup();
    
    let sketch: Vec<u64> = kmer_hashes
        .into_iter()
        .take(sketch_size)
        .collect();

    if sketch.is_empty() {
        anyhow::bail!("Failed to create sketch");
    }

    Ok(sketch)
}

/// Calculates reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
        .collect()
}

/// Converts Jaccard index to Mash distance
fn jaccard_to_mash(jaccard: f64, k: usize) -> f64 {
    if jaccard <= 0.0 {
        return 1.0;
    }
    if jaccard >= 1.0 {
        return 0.0;
    }
    
    // Mash distance formula: d = -1/k * ln(2*J/(1+J))
    let numerator = 2.0 * jaccard;
    let denominator = 1.0 + jaccard;
    let ratio = numerator / denominator;
    
    if ratio <= 0.0 {
        return 1.0;
    }
    
    let distance = -1.0 / (k as f64) * ratio.ln();
    distance.max(0.0).min(1.0)
}

/// Calculates a signature (hash) for a sketch based on its MinHash values
fn calculate_sketch_signature(sketch_values: &[u64]) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let mut hasher = DefaultHasher::new();
    sketch_values.hash(&mut hasher);
    hasher.finish()
}

/// Manages multiple sketches indexed by signature
pub struct SketchManager {
    pub sketches: HashMap<u64, Sketch>,      // Índice principal por signature
    pub name_to_sig: HashMap<String, u64>,   // Índice auxiliar name → signature
    pub kmer_size: usize,
    pub sketch_size: usize,
}

impl SketchManager {
    /// Creates a new empty SketchManager
    pub fn new(kmer_size: usize, sketch_size: usize) -> Self {
        SketchManager {
            sketches: HashMap::new(),
            name_to_sig: HashMap::new(),
            kmer_size,
            sketch_size,
        }
    }

    /// Adds a sketch to the manager
    pub fn add_sketch(&mut self, sketch: Sketch) -> Result<()> {
        let signature = sketch.signature;
        let name = sketch.name.clone();
        
        self.sketches.insert(signature, sketch);
        self.name_to_sig.insert(name, signature);
        
        Ok(())
    }

    // ═══════════════════════════════════════════════════════════════
    // FUNCIONES PRINCIPALES (usan signature)
    // ═══════════════════════════════════════════════════════════════

    /// Gets a sketch by signature
    pub fn get_sketch(&self, signature: u64) -> Option<&Sketch> {
        self.sketches.get(&signature)
    }

    /// Checks if a signature exists
    pub fn contains(&self, signature: u64) -> bool {
        self.sketches.contains_key(&signature)
    }

    /// Gets all signatures
    pub fn get_all_signatures(&self) -> Vec<u64> {
        self.sketches.keys().copied().collect()
    }

    // ═══════════════════════════════════════════════════════════════
    // FUNCIONES AUXILIARES (usan name, con sufijo _by_name)
    // ═══════════════════════════════════════════════════════════════

    /// Gets a sketch by name
    pub fn get_sketch_by_name(&self, name: &str) -> Option<&Sketch> {
        self.name_to_sig.get(name)
            .and_then(|sig| self.sketches.get(sig))
    }

    /// Checks if a name exists
    pub fn contains_by_name(&self, name: &str) -> bool {
        self.name_to_sig.contains_key(name)
    }

    /// Gets signature by name
    pub fn get_signature_by_name(&self, name: &str) -> Option<u64> {
        self.name_to_sig.get(name).copied()
    }

    // ═══════════════════════════════════════════════════════════════
    // FUNCIONES GENERALES
    // ═══════════════════════════════════════════════════════════════

    /// Returns the number of sketches
    pub fn length(&self) -> usize {
        self.sketches.len()
    }

    /// Loads sketches from a directory of FASTA files
    pub fn load_from_directory(&mut self, dir_path: &Path) -> Result<usize> {
        let mut count = 0;

        for entry in std::fs::read_dir(dir_path)
            .with_context(|| format!("Failed to read directory: {}", dir_path.display()))?
        {
            let entry = entry.context("Failed to read directory entry")?;
            let path = entry.path();

            if path.is_file() {
                let extension = path.extension().and_then(|s| s.to_str());
                if matches!(extension, Some("fasta") | Some("fa") | Some("fna")) {
                    match Sketch::new(&path, self.kmer_size, self.sketch_size) {
                        Ok(sketch) => {
                            self.add_sketch(sketch)?;
                            count += 1;
                        }
                        Err(e) => {
                            eprintln!("Warning: Failed to create sketch from {}: {}", 
                                     path.display(), e);
                        }
                    }
                }
            }
        }

        Ok(count)
    }

    /// Calculates pairwise distances for all sketches
    pub fn pairwise_distances(&self) -> Vec<(u64, u64, f64)> {
        let mut results = Vec::new();
        let signatures: Vec<u64> = self.sketches.keys().copied().collect();

        for i in 0..signatures.len() {
            for j in (i + 1)..signatures.len() {
                let sig1 = signatures[i];
                let sig2 = signatures[j];
                
                if let (Some(sketch1), Some(sketch2)) = 
                    (self.sketches.get(&sig1), self.sketches.get(&sig2)) {
                    let distance = sketch1.distance(sketch2);
                    results.push((sig1, sig2, distance));
                }
            }
        }

        results
    }
}

/// Calculates distances from one query against many reference sketches
pub fn pw_one_to_many(
    query_manager: &SketchManager,
    ref_manager: &SketchManager,
    subset_signatures: &[u64],
    threshold: f64,
) -> Vec<(u64, u64, f64)> {
    let mut results = Vec::new();

    // Obtener el único sketch del query_manager
    let query_sketch = if let Some(sketch) = query_manager.sketches.values().next() {
        sketch
    } else {
        return results;
    };

    let query_sig = query_sketch.signature;

    // Iterar sobre los signatures de referencia
    for &ref_sig in subset_signatures {
        if let Some(ref_sketch) = ref_manager.get_sketch(ref_sig) {
            let distance = query_sketch.distance(ref_sketch);
            
            if distance >= threshold {
                results.push((query_sig, ref_sig, distance));
            }
        }
    }

    results
}

/// Loads SketchManager from database
pub fn load_sketch_manager_from_db(
    conn: &duckdb::Connection,
    kmer_size: usize,
    sketch_size: usize,
) -> anyhow::Result<SketchManager> {
    use crate::db::db::get_all_sketch_objects;

    let mut manager = SketchManager::new(kmer_size, sketch_size);
    let sketches = get_all_sketch_objects(conn)
        .context("Error loading sketches from database")?;

    for (_sample_name, sketch) in sketches {
        manager.add_sketch(sketch)
            .context("Error adding sketch to manager")?;
    }

    Ok(manager)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("GCTA"), "TAGC");
    }

    #[test]
    fn test_sketch_creation() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, ">test_sequence").unwrap();
        writeln!(temp_file, "ATCGATCGATCGATCG").unwrap();

        let sketch = Sketch::new(temp_file.path(), 3, 10).unwrap();
        assert_eq!(sketch.kmer_size, 3);
        assert_eq!(sketch.sketch_size, 10);
        assert!(sketch.signature != 0);
        assert!(!sketch.sketch.is_empty());
    }

    #[test]
    fn test_sketch_manager() {
        let mut manager = SketchManager::new(3, 10);

        let mut temp_file1 = NamedTempFile::new().unwrap();
        writeln!(temp_file1, ">seq1").unwrap();
        writeln!(temp_file1, "ATCGATCGATCGATCG").unwrap();

        let sketch1 = Sketch::new(temp_file1.path(), 3, 10).unwrap();
        let sig1 = sketch1.signature;
        let name1 = sketch1.name.clone();
        
        manager.add_sketch(sketch1).unwrap();

        assert_eq!(manager.length(), 1);
        assert!(manager.contains(sig1));
        assert!(manager.get_sketch(sig1).is_some());
        assert!(manager.contains_by_name(&name1));
        assert!(manager.get_sketch_by_name(&name1).is_some());
        assert_eq!(manager.get_signature_by_name(&name1), Some(sig1));
    }

    #[test]
    fn test_distance_calculation() {
        let mut temp_file1 = NamedTempFile::new().unwrap();
        writeln!(temp_file1, ">seq1").unwrap();
        writeln!(temp_file1, "ATCGATCGATCGATCG").unwrap();

        let mut temp_file2 = NamedTempFile::new().unwrap();
        writeln!(temp_file2, ">seq2").unwrap();
        writeln!(temp_file2, "ATCGATCGATCGATCG").unwrap();

        let sketch1 = Sketch::new(temp_file1.path(), 3, 10).unwrap();
        let sketch2 = Sketch::new(temp_file2.path(), 3, 10).unwrap();

        let distance = sketch1.distance(&sketch2);
        assert!(distance > 0.9); // Should be very similar (same sequence)
    }

    #[test]
    fn test_pw_one_to_many() {
        let mut query_manager = SketchManager::new(3, 10);
        let mut ref_manager = SketchManager::new(3, 10);

        let mut temp_query = NamedTempFile::new().unwrap();
        writeln!(temp_query, ">query").unwrap();
        writeln!(temp_query, "ATCGATCGATCGATCG").unwrap();

        let mut temp_ref1 = NamedTempFile::new().unwrap();
        writeln!(temp_ref1, ">ref1").unwrap();
        writeln!(temp_ref1, "ATCGATCGATCGATCG").unwrap();

        let query_sketch = Sketch::new(temp_query.path(), 3, 10).unwrap();
        let ref_sketch1 = Sketch::new(temp_ref1.path(), 3, 10).unwrap();
        
        let ref_sig1 = ref_sketch1.signature;

        query_manager.add_sketch(query_sketch).unwrap();
        ref_manager.add_sketch(ref_sketch1).unwrap();

        let subset_sigs = vec![ref_sig1];
        let distances = pw_one_to_many(&query_manager, &ref_manager, &subset_sigs, 0.0);

        assert_eq!(distances.len(), 1);
        assert!(distances[0].2 > 0.9);
    }

    #[test]
    fn test_jaccard_to_mash() {
        // Test casos extremos
        assert_eq!(jaccard_to_mash(0.0, 21), 1.0);
        assert_eq!(jaccard_to_mash(1.0, 21), 0.0);
        
        // Test caso intermedio
        let distance = jaccard_to_mash(0.5, 21);
        assert!(distance > 0.0 && distance < 1.0);
    }
}
