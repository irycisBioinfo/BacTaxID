# BacTaxID - Fast Bacterial Taxonomy Identification Tool

<div align="center">

  <img src="images/BacTaxID_Logo.png" alt="BacTaxID Logo" width="400"/>

  <br/>

  <em>High-performance bacterial sub-genus classification using advanced sketching algorithms</em>

</div>

## Overview

**BacTaxID** is a universal, genome-based bacterial typing system that overcomes the limitations of traditional species-specific frameworks by unifying scalable classification with biologically meaningful nomenclature. Developed in Rust, it uses Binwise Densified MinHash combined with ntHash for rapid k-mer sketching, achieving highly efficient $O(S)$ distance comparison complexity while maintaining strict proportionality to Average Nucleotide Identity (ANI).

Applied to **2.3 million genomes** across 67 genera from the *All the Bacteria* database, BacTaxID demonstrates:

* Universal concordance with species and sub-species classification systems.


* Epidemiologically relevant resolution at intermediate tiers (L₃ ≈ 99.5% ANI ≈ MLST).


* Sub-clonal outbreak detection capabilities at the finest tier (L₅ ≈ 99.99% ANI ≈ 3-29 SNPs/Mb).


* Scalability to millions of genomes through sublinear hierarchical search complexity ($O(\log N)$ to $O(N^{0.5})$).



---

## Installation & Building from Source

BacTaxID is implemented in Rust, utilizing the `rayon` library for efficient multi-threaded parallelization.

```bash
# Clone the repository
git clone https://github.com/irycisBioinfo/BacTaxID.git
cd BacTaxID

# Ensure your Rust toolchain is up to date
rustup update

# Build the release binary
cargo build --release

# The compiled binary will be available at:
./target/release/bactaxid --help

```

### Troubleshooting Build Issues

If you encounter an error parsing `Cargo.lock` during compilation (e.g., `error: failed to parse lock file at .../Cargo.lock`), this is typically due to an outdated local Rust toolchain or a lockfile conflict in shared HPC environments.

To resolve this, remove the existing lockfile and re-compile to regenerate it cleanly:

```bash
rm Cargo.lock
cargo build --release

```

---

## Quick Start & Commands

BacTaxID operates via four primary subcommands, encapsulating its multi-tiered architecture within a portable, single-file DuckDB database. Configuration parameters are passed directly via CLI flags.

### 1. `bactaxid init`

Initializes an empty portable DuckDB database with the specified hierarchical parameters.

```bash
# Initialize with standard defaults (6 levels, k=31, sketch=3000)
./target/release/bactaxid init --db bacteria.db

# Or initialize with custom epidemiological parameters
./target/release/bactaxid init \
  --db bacteria.db \
  --kmer-size 31 \
  --sketch-size 3000 \
  --click-size 5 \
  --click-threshold 0.8 \
  --reference-size 100

```

### 2. `bactaxid update`

Populates or dynamically expands an existing database by sketching assemblies and performing iterative pseudo-clique cluster assignments.

```bash
./target/release/bactaxid update \
  --db bacteria.db \
  --fasta-dir /path/to/genomes/ \
  --threads 16

```

* Generates sketches and computes deterministic signatures.


* Determines **Classifier** (core reference anchors) or **Satellite** (peripheral members) status to prevent single-linkage chaining artifacts.


* Seamlessly incorporates novel diversity without disrupting established historical classifications.



### 3. `bactaxid classify`

Queries new or unknown FASTA assemblies against an established database offline, assigning hierarchical typing codes deterministically from L₀ down to L₅.

```bash
./target/release/bactaxid classify \
  --db bacteria.db \
  --query query_genome.fasta

```

### 4. `bactaxid distance`

Calculates sketch Jaccard distances and transformed whole-genome ANI between assemblies.

```bash
./target/release/bactaxid distance \
  --genome1 genomeA.fasta \
  --genome2 genomeB.fasta \
  --kmer-size 31 \
  --sketch-size 3000

```

---

## Technical Architecture & Algorithm Details

### High-Performance Sketching Engine

* **ntHash Streaming**: Generates canonical k-mer hashes on-the-fly via a rolling window, avoiding memory-intensive k-mer materialization.


* **Binwise Densified MinHash**: Partitions the 64-bit hash space into predefined bins, retaining only the minimum hash value per bin in constant $O(1)$ time.


* **Distance Computation**: Exploits direct element-wise iteration across aligned binwise arrays, comparing minimum hash values at matching bin indices in $O(S)$ time per pairwise comparison (without the overhead of sort-traversal algorithms).


* **ANI Conversion**: The relationship between Jaccard similarity and ANI is established through the Mash distance formula, maintaining strong linearity for ANI values ≥85%.



### Hierarchical Classification Algorithm

1. **Top-Down Hierarchical Pruning:** Queries are evaluated level-by-level (L₀ to L₅). At each tier, the search space is restricted exclusively to core classifiers located within the matching parent container, instantly pruning disjoint branches.


2. **Representative Node Filtering:** Distances are evaluated exclusively against **Classifying Nodes** (the core seeds of the pseudo-clique). Peripheral **Satellite Nodes** are excluded from serving as references to prevent chaining artifacts.


3. **De Novo Group Formation:** When no existing cluster matches, the system initiates graph-theoretic clique detection to identify complete subgraphs (maximal cliques) where all members satisfy pairwise distance thresholds, establishing highly cohesive novel clusters.


4. **Explicit Singleton State:** Isolates failing similarity or density thresholds retain an explicit unassigned state ('0'), avoiding forced classification of biologically unsupported singletons.



---

## Pre-computed Resources & Database Integration

To establish an initial global standard, we have pre-computed schemes spanning 67 genera and >2.3 million genomes.

* **Web Platform**: Interactive exploration tools and KronaPlots are freely available at [www.bactaxid.org](http://www.bactaxid.org).


* **Data Repository**: Full classification datasets are deposited at [Zenodo](https://zenodo.org/records/17791772).


* **Metadata Traceability**: All BacTaxID databases share identical unique primary keys (NCBI/EBI BioSample accessions) with the *All the Bacteria* resource, enabling seamless SQL integration with comprehensive antimicrobial, biocide, and metal resistance gene profiles.



---

## Citation

If you use BacTaxID in your research, please cite:

**Díez Fernández de Bobadilla M, Fernández Lanza V (2025). BacTaxID: A universal framework for standardized bacterial classification.**

## License

GPL-3.0 • **Maintainers**: Val F. Lanza, Miguel Diez Fernandez de Bobadilla

**Issues**: [GitHub Issues](https://github.com/irycisBioinfo/BacTaxID/issues)
