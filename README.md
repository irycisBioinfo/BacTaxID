<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# BacTaxID - Fast Bacterial Taxonomy Identification Tool

## Overview

**BacTaxID** is a high-performance bacterial sub-genus classification tool developed in Rust that uses advanced sketching algorithms for rapid and accurate classification of bacterial genomes. The tool implements a BinDash-like algorithm, a state-of-the-art sketching technique derived from BinHash, to provide efficient similarity with superior speed and accuracy compared to traditional methods.

## Key Features

### 🚀 **High-Performance Sketching**

- **BinDash-like Algorithm**: O(log S) Complexity: Uses direct binning instead of heap operations for optimal performance
- **Rolling-Hash Streaming**: Efficient processing of FASTA files with ntHash (https://github.com/bcgsc/ntHash)
- **Parallel Processing**: Leverages Rayon for multi-threaded sketch comparison and analysis


### 🗄️ **Optimized Database Architecture**
- **DuckDB Backend**: Fast analytical database optimized for genomic data queries
- **Dynamic Level Support**: Configurable taxonomic levels with custom ANI thresholds
- **Sketch Serialization**: Efficient storage and retrieval of sketch objects using bincode


### 🧬 **Advanced Taxonomic Classification**

- **Hierarchical Classification**: Multi-level taxonomic assignment (L_0 to L_N)
- **ANI-based Thresholds**: Average Nucleotide Identity calculations for similarity assessment
- **Graph-based Analysis**: Clique detection using PetGraph for community identification


### 📊 **Intelligent Clustering**

- **Clique Detection**: Identifies maximal cliques for taxonomic grouping
- **Connectivity Analysis**: Distance-based edge filtering with configurable thresholds
- **Community Detection**: Automatic discovery of taxonomic clusters
- **Reference Management**: Dynamic reference database with size-based optimization


## Technical Architecture

### Core Components

1. **Sketching Engine** (`src/sketch/sketching.rs`)
    - BinDash sketch generation from FASTA files
    - Jaccard distance calculation converted to ANI
    - SketchManager for batch operations
2. **Database Layer** (`src/db/db.rs`)
    - Unified taxonomy storage schema
    - Batch operations for high-throughput analysis
    - Automatic indexing and optimization
3. **Graph Analysis** (`src/graph/graph.rs`)
    - PetGraph-based clique detection
    - Edge management and filtering
    - Community analysis algorithms
4. **Command Interface** (`src/commands/`)
    - CLI-based workflow management
    - TOML configuration support
    - Batch processing capabilities

### Configuration

BacTaxID uses TOML configuration files for flexible parameter management:

```toml
genus = "Escherichia"
acronym = "EC"
levels = "[0.95, 0.98, 0.99, 0.999, 0.9999]"
kmer_size = 21
sketch_size = 1000
click_size = 8
click_threshold = 0.8
reference_size = 100
```

## Requirements

- **Rust**: 2021 edition or later
- **System Memory**: Minimum 2GB RAM (configurable)
- **Storage**: SSD recommended for optimal I/O performance
- **CPU**: Multi-core processor recommended for parallel processing


## Installation

```bash
git clone https://github.com/irycisBioinfo/BacTaxID
cd BacTaxID
cargo build --release
```


## Quick Start

1. **Initialize Database**:

```bash
./target/release/bactaxid init --config config.toml --db bacteria.db
```

2. **Add Reference Genomes**:

```bash
./target/release/bactaxid update --db bacteria.db --files reference_list.txt
```

3. **Classify New Samples**:

```bash
./target/release/bactaxid classify --db bacteria.db --query sample.fasta
```


## Algorithm Details

### BinDash Sketching

BacTaxID implements the BinDash-like algorithm, which provides several advantages over traditional MinHash approaches:

- **Direct Binning**: Eliminates heap operations for O(log S) complexity
- **Better Accuracy**: Improved estimation of Jaccard similarities
- **Memory Efficiency**: Reduced memory footprint compared to traditional sketching


### Taxonomic Assignment

The tool uses a hierarchical approach with configurable ANI thresholds:

1. **Hierarchical Classification**: Build self-explanatory classification in N ANI levels.
2. **Best Hit Analysis**: Identifies closest matches in reference database.
3. **Robust Reference DataBase**: Identifies reference communities avoiding chain-effect.
4. **Confidence Assessment**: Evaluates classification reliability.
5. **Clique Detection**: Groups similar organisms into taxonomic clusters.



## Contributing

We welcome contributions to BacTaxID! Please see our contributing guidelines and feel free to submit issues or pull requests.

## License

BacTaxID is released under the GPL-3.0 license. See `LICENSE` file for details.

## Citation

If you use BacTaxID in your research, please cite:

*[Citation information will be added upon publication]*

## Contact

- **Development Team**: irycisBioinfo
- **Maintainers**: Val F. Lanza, Miguel Diez Fernandez de Bobadilla
- **Issues**: Please use GitHub Issues for bug reports and feature requests

***

**BacTaxID** - Empowering bacterial taxonomy identification through advanced computational methods and optimized performance.

