# BacTaxID - Fast Bacterial Taxonomy Identification Tool

<div align="center">
  <img src="images/BacTaxID_Logo.png" alt="BacTaxID Logo" width="400"/>
  <br/>
  <em>High-performance bacterial sub-genus classification using advanced sketching algorithms</em>
</div>

## Overview

**BacTaxID** is a universal, genome-based bacterial typing system that overcomes limitations of traditional species-specific frameworks. Developed in Rust, it uses **Binwise Densified MinHash** (BinDash) combined with **ntHash** for rapid k-mer sketching, achieving O(log S) complexity while maintaining strict proportionality to Average Nucleotide Identity (ANI).

Applied to **2.3 million genomes** across 67 genera, BacTaxID demonstrates:
- Universal concordance with species/sub-species classifications
- Epidemiologically relevant resolution (L₃ ≈ 99% ANI ≈ MLST)
- Sub-clonal outbreak detection (L₅ ≈ 99.99% ANI ≈ 5-27 SNPs/Mb)
- Scalability to millions of genomes with hierarchical search complexity

## Key Features

### 🚀 **High-Performance Sketching**
- **BinDash Algorithm**: Direct binning eliminates heap operations (O(log S) vs O(S log S))
- **ntHash Streaming**: Rolling-hash k-mer decomposition without materialization
- **Parallel Processing**: Rayon-powered multi-threaded distance computations

### 🗄️ **Optimized Database Architecture**
- **Signature-based Indexing**: xxHash64 u64 signatures for unique genome identification
- **UBIGINT[] Sketch Storage**: Native SQL arrays to universal usage
- **Automatic Duplicate Detection**: `duplicates` table links identical sketches
- **DuckDB Backend**: Self-contained portable database with SQL/Python/R/CLI APIs

### 🧬 **Advanced Typing System**
- **Pseudo-Clique Clustering**: Prevents chaining artifacts via click_threshold criterion
- **Hierarchical Resolution**: L₀ (96% ANI) → L₅ (99.99% ANI) configurable levels
- **Classifier/Satellite States**: Robust reference management avoiding chain effects
- **Monotonicity Constraints**: Ensures hierarchical consistency across levels

### 📊 **Intelligent Clustering**
- **Maximal Clique Detection**: PetGraph-based complete subgraph identification
- **ANI-based Thresholds**: Direct quantitative link to Average Nucleotide Identity
- **Reference Community Detection**: Size-based optimization (reference_size parameter)
- **Dynamic Scheme Evolution**: Incremental updates without full re-clustering

## Technical Architecture

### Core Components

1. **Sketching Engine** (`src/sketch/sketching.rs`)
   - BinDash sketch generation: partition hash space into bins, retain minimum per bin
   - Jaccard similarity → ANI conversion: `ANI = 1 - (-(1/k)*ln(2J/(1+J)))`
   - xxHash64 signature generation for unique genome identification

2. **Database Layer** (`src/db/db.rs`)
   - `sketches` table: signature (UBIGINT PRIMARY KEY), UBIGINT[] sketch array
   - `duplicates` table: maps duplicate samples to canonical signatures
   - `code` table: hierarchical classifications with L_i_int, L_i_full, L_i_state columns

3. **Graph Analysis** (`src/graph/graph.rs`)
   - Maximal clique detection using PetGraph's `maximal_cliques` algorithm
   - Edge filtering by ANI threshold and parent-level group consistency
   - Strategic edge pruning post-clique formation for next-level searches

4. **Command Interface** (`src/commands/`)
   - `init`: Initialize database with TOML configuration
   - `update`: Add reference genomes, detect duplicates, build hierarchy
   - `classify`: Read-only classification of query genomes
   - `distance`: Pairwise ANI distance matrix generation

### Configuration (`config.toml`)

```toml
genus = "Escherichia"
acronym = "EC"
levels = "[0.96, 0.98, 0.99, 0.995, 0.999, 0.9999]"  # ANI thresholds
kmer_size = 31
sketch_size = 4000
click_size = 5         # Minimum clique size for new cluster
click_threshold = 0.8  # Fraction of cluster members query must match
reference_size = 100   # Maximum classifiers per cluster
```

## Commands

### `bactaxid init`
Initialize database with taxonomic levels and metadata.

```bash
bactaxid init --config config.toml --db bacteria.db
```

### `bactaxid update`
```bash
bactaxid update --db bacteria.db --files references.txt [--cpus 8] [--debug]
```
- Processes FASTA files sequentially, generates sketches, computes signatures
- Detects duplicates: identical signatures → `duplicates` table entry
- Hierarchical classification: best-hit search → clique detection if no assignment
- **Classifier/Satellite assignment**: based on click_threshold and reference_size
- **Debug mode**: saves all pairwise distances to `debug` table

### `bactaxid classify` 
```bash
bactaxid classify --db bacteria.db --queries queries.txt --output results.tsv [--verbose]
```
- **Read-only** classification against pre-built reference scheme
- **Hierarchical search**: compares query only to classifiers sharing parent-level code
- **Early stopping**: terminates at first level without valid best-hit
- **Output**: TSV with query_id, signature, best_hit, similarity, levels_reached, taxonomic codes

### `bactaxid distance` 
```bash
bactaxid distance --db bacteria.db --ids samples.txt --output distances.phylip \
  --format phylip [--threshold 0.95] [--threads 8]
```
- **Pairwise ANI distances**: computes all-vs-all similarity matrix
- **Formats**:
  - `tsv`/`csv`: Long format (sample1, sample2, ANI, distance)
  - `phylip`: Square matrix for phylogenetic analysis
- **Duplicate handling**: transparently resolves samples to canonical signatures
- **Parallel computation**: Rayon-powered multi-threaded distance calculations
- **Threshold filtering** (TSV/CSV only): reports only pairs ≥ threshold

## Algorithm Details

### BinDash Sketching

1. **Hash Space Partitioning**: Divide 2⁶⁴ hash space into `sketch_size` bins
2. **Streaming k-mer Processing**: ntHash generates canonical k-mer hashes on-the-fly
3. **Bin Assignment**: `bin_index = hash / bin_size`, retain minimum hash per bin
4. **Sketch Comparison**: Count matching bin values → Jaccard estimate `J = matches / S`
5. **ANI Conversion**: `ANI = 1 - (-(1/k)*ln(2J/(1+J)))`, accurate for ANI ≥ 85%

**Advantages**:
- **O(1) per k-mer** (vs O(log S) for bottom-k MinHash heap operations)
- **Memory-efficient**: Fixed-size vector regardless of genome length
- **ANI proportionality**: Direct quantitative link to evolutionary distance

### Hierarchical Classification Algorithm

#### **Phase 1: Best-Hit Search (Classifiers)**
```
For level i:
  1. Retrieve classifiers C_i sharing parent code (level i-1)
  2. Compute distances: query → all c ∈ C_i
  3. If best_hit ∈ C_i satisfies threshold AND
     click_criterion(query, C_i, best_hit.code):
       → Assign query to best_hit.code
       → Promote to Classifier if reference_size not exceeded
```

**Click Criterion** (prevents chaining):
```python
def click_criterion(query, classifiers, code):
    cluster = [c for c in classifiers if c.code == code]
    matches = sum(1 for c in cluster if distance(query, c) >= threshold)
    return matches / len(cluster) >= click_threshold
```

#### **Phase 2: Clique Detection (De Novo Clustering)**
```
If no valid best-hit:
  1. Retrieve ALL genomes G sharing parent code (classifiers + satellites)
  2. Construct distance graph: nodes = G, edges = pairs meeting threshold
  3. Find maximal cliques using PetGraph::maximal_cliques()
  4. If clique size >= click_size:
       → Create new code (next integer for parent group)
       → Assign ALL clique members as Classifiers
       → Assign non-clique graph members as Satellites
       → Prune edges for next-level search
```

**Monotonicity Constraint**: `code_full[i]` must start with `code_full[i-1]`

### Signature System (xxHash64)

- **Unique Identification**: Each sketch assigned 64-bit signature based on content hash
- **Duplicate Detection**: Identical sketches share signature → `duplicates` table
- **Primary Key**: `sketches.signature` replaces `sketches.sample` for scalability
- **Metadata Tracking**: `duplicates` maps sample names → canonical signatures

## Recent Updates (Dec 2025)

| Feature | Description | Benefit |
|---------|-------------|---------|
| **Signature Architecture** | xxHash64 u64 identifiers, `duplicates` table | Efficient duplicate handling, robust indexing |
| **UBIGINT[] Storage** | Native SQL sketch arrays | Direct SQL queries, no deserialization |
| **`classify` Command** | Read-only hierarchical classification | Fast prospective typing, no DB modifications |
| **`distance` Command** | Pairwise ANI matrices (TSV/CSV/PHYLIP) | Phylogenetic analysis, outbreak investigation |
| **Parallel Distance** | Rayon multi-threading | Minutes for millions of comparisons |
| **English Documentation** | Full code comments translation | Improved accessibility |

## Quick Start

```bash
# 1. Clone and build
git clone https://github.com/irycisBioinfo/BacTaxID
cd BacTaxID
cargo build --release

# 2. Initialize database
./target/release/bactaxid init --config config.toml --db escherichia.db

# 3. Build reference scheme
./target/release/bactaxid update --db escherichia.db --files references.txt --cpus 8

# 4. Classify new samples
./target/release/bactaxid classify --db escherichia.db --queries queries.txt \
  --output classifications.tsv --verbose

# 5. Generate distance matrix
./target/release/bactaxid distance --db escherichia.db --ids samples.txt \
  --output distances.phylip --format phylip --threads 8
```

## Performance Characteristics

- **Sketch generation**: ~1-2s per genome (5 Mb, k=31, S=4000)
- **Classification**: O(log N) per level due to hierarchical search
- **All-vs-all distances**: O(N²S) parallelized with Rayon
- **Memory**: ~2GB minimum, scales with database size
- **Database size**: ~50MB per 1000 genomes (k=31, S=4000)

## Validation (2.3M Genomes, 67 Genera)

- **Species concordance**: L₀-L₁ (96-98% ANI) ≈ 95% ANI species threshold
- **MLST equivalence**: L₃ (99% ANI) universal peak across 11 genera
- **Outbreak detection**: L₅ (99.99% ANI) → 5.92-27 SNPs/Mb (genus-dependent)
- **Completeness**: 98% at L₀, variable at L₅ (population structure-dependent)

## Requirements

- **Rust**: 2021 edition or later
- **System Memory**: 2GB minimum (scales with dataset)
- **Storage**: SSD recommended for I/O performance
- **CPU**: Multi-core beneficial for parallel operations

## Installation

```bash
git clone https://github.com/irycisBioinfo/BacTaxID
cd BacTaxID
cargo build --release
```

**Pre-compiled binaries** and **67 pre-computed genus schemes** available at:
- 🌐 **www.bactaxid.org** (interactive tools, KronaPlots)
- 📦 **zenodo.org/doi/10.5281/zenodo.17226627** (classification data)

## Citation

If you use BacTaxID in your research, please cite:

**Díez Fernández de Bobadilla M, Fernández Lanza V (2025). BacTaxID: A universal framework for standardized bacterial classification. _[Manuscript in preparation]_**

## License

GPL-3.0 • **Maintainers**: Val F. Lanza, Miguel Diez Fernandez de Bobadilla  
**Issues**: [GitHub Issues](https://github.com/irycisBioinfo/BacTaxID/issues)

---

**BacTaxID** - Empowering bacterial taxonomy identification through advanced computational methods, hierarchical resolution, and universal standardization.