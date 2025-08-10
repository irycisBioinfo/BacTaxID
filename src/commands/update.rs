/// Arguments for the update subcommand
#[derive(Args, Debug)]
pub struct UpdateArgs {
    /// Path to the DuckDB database
    #[arg(long, required = true, value_name = "DB_PATH")]
    pub db: String,

    /// Number of CPUs to parallelize with rayon
    #[arg(long, value_name = "N_CPUS", default_value_t = 1)]
    pub cpus: usize,

    /// Plain text file with paths to FASTA files (one per line)
    #[arg(long, required = true, value_name = "FILES_LIST")]
    pub files: String,

    #[arg(long, required= false, value_name = "DEBUG")]
    pub debug: bool
}

/// Context with preloaded information needed by the entire `update` flow
pub struct UpdateCtx<'a> {
    /// Mutable connection to DuckDB (reusable in all phases)
    pub conn: &'a mut Connection,
    /// Levels table (L_0 … L_N → dist) as an ordered vector
    pub levels: Vec<(String, f64)>,
    /// Entire metadata row as (column → value in UTF-8)
    pub metadata_map: HashMap<String, String>,
}

impl<'a> UpdateCtx<'a> {
    /// Convenience accessors for metadata
    pub fn kmer_size(&self) -> usize {
        self.metadata_map["kmer_size"].parse().unwrap_or(21)
    }

    pub fn sketch_size(&self) -> usize {
        self.metadata_map["sketch_size"].parse().unwrap_or(1000)
    }

    pub fn click_size(&self) -> usize {
        self.metadata_map["click_size"].parse().unwrap_or(50)
    }

    pub fn click_threshold(&self) -> f64 {
        self.metadata_map["click_threshold"].parse().unwrap_or(0.025)
    }

    pub fn genus(&self) -> &str {
        &self.metadata_map["genus"]
    }

    pub fn acronym(&self) -> &str {
        &self.metadata_map["acronym"]
    }

    pub fn levels_str(&self) -> &str {
        &self.metadata_map["levels"]
    }

    /// Access to the mutable connection
    pub fn connection_mut(&mut self) -> &mut Connection {
        self.conn
    }

    pub fn reference_size(&self) -> usize {
        self.metadata_map["reference_size"].parse().unwrap_or(100)
    }

    /// Number of levels
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    /// Get distance for a specific level
    pub fn get_level_distance(&self, level: &str) -> Option<f64> {
        self.levels.iter().find(|(l, _)| l == level).map(|(_, d)| *d)
    }
}

/// Reads the `levels` table and returns a Vec<(String, f64)> ordered by level
fn load_levels_vec(conn: &Connection) -> Result<Vec<(String, f64)>> {
    // ...
}

/// Reads metadata (a single row) → HashMap<column, value_as_text>
fn load_metadata_map(conn: &Connection) -> Result<HashMap<String, String>> {
    // ...
}

/// Main function for the update command
pub fn update_command(args: &UpdateArgs) -> Result<()> {
    // ...
    // TODO: Here you will add logic for:
    // 1. Compare query vs reference sketches
    // 2. Update counters in level_counts
    // 3. Insert/update records in code, code_full, code_state tables
    // 4. Save updated SketchManager
    // ...
}

pub struct Query {
    pub sample_name: String,
    pub code: Vec<usize>,         // vector of positive integers
    pub code_full: Vec<String>,   // vector of strings
    pub code_state: Vec<String>,  // vector of strings
    pub sketch: SketchManager,
}

impl Query {
    /// Creates a new Query loading the SketchManager+Sketch, and vectors of the appropriate size
    ///
    /// # Arguments
    /// * `path` - Path to the input FASTA file.
    /// * `ctx` - Context that must be able to (for example) know the number of levels.
    pub fn new(path: &Path, ctx: &UpdateCtx) -> anyhow::Result<Self> {
        // ...
    }
}

/// Processes **one** FASTA file:
/// 1. Creates a `Sketch` of the sample using `Sketch::new` (which reads the FASTA internally).
/// 2. Creates the `query_class` vector with the length of the levels.
/// 3. Returns `(sample_name, query_sketch, query_class)`
///
/// * `query_sketch` is the sketch created from the FASTA file.
/// * `query_class` is a vector of `u64` with length `ctx.levels_map.len()`
///   initialized to 0; it will be used later to count matches per level.
pub fn update_single_file(
    fasta_path: &Path,
    ctx: &mut UpdateCtx,
    _sketch_manager: &mut SketchManager,
    debug: bool
) -> Result<()> {
    // ...
}

/// Retrieves classifiers that meet specific level and group conditions.
/// 
/// # Arguments
/// * `conn` - Connection to the DuckDB database
/// * `level` - Taxonomic level (0, 1, 2, 3, 4...)
/// * `group` - Taxonomic group to filter (ignored if level == 0)
/// 
/// # Returns
/// Vector of tuples (sample, code, code_full, code_state) where:
/// - code_state.L_{level} == "C" 
/// - code_full.L_{level} == group (except if level == 0)
/// 
/// # Example
/// ```
/// // Get all classifiers at level 0 (ignores group)
/// let results = retrieve_classifiers(conn, 0, "0")?;
/// 
/// // Get classifiers at level 2 for group "Escherichia"
/// let results = retrieve_classifiers(conn, 2, "1.1")?;
/// ```
pub fn retrieve_classifiers(
    conn: &Connection,
    level: usize,
    group: &str,
    condition: &str,
) -> Result<Vec<(String, usize, String, String)>> {
    // ...
}

/// Finds the maximum similarity and returns the corresponding tuple from ref_db.
/// 
/// # Arguments
/// * `distances` - Vector of tuples (query_name, ref_name, distance)
/// * `ref_db` - Vector of tuples (sample, code, code_full, code_state) from retrieve_classifiers
/// 
/// # Returns
/// Option with the tuple from ref_db that has the sample corresponding to the best hit,
/// or None if no match is found or there are no valid distances.
pub fn best_hit(
    distances: &[(String, String, f64)],
    ref_db: &[(String, usize, String, String)],
) -> Option<(String, usize, String, String)> {
    // ...
}

/// Determines if a best_hit passes the classification criterion.
///
/// # Arguments
/// * `distances`       - Vector of tuples (query_name, ref_name, dist).
/// * `best_hit`        - Tuple (sample, code, code_full, code_state) of the best hit.
/// * `ref_db`          - Vector of tuples (sample, code, code_full, code_state).
/// * `click_threshold` - Minimum proportion threshold (e.g. 0.025).
/// * `reference_size`  - Maximum allowed group size; if matched or exceeded, fails.
///
/// # Returns
/// `true` if (count_dist / count_ref) >= click_threshold **and** count_ref < reference_size;
/// `false` otherwise.
pub fn is_classifier(
    distances: &[(String, String, f64)],
    best_hit: &(String, usize, String, String),
    ref_db: &[(String, usize, String, String)],
    click_threshold: f64,
    reference_size: usize,
) -> bool {
    // ...
}

/// Searches for cliques in specific taxonomic levels and updates the Query object.
///
/// # Arguments
/// * `conn` - Mutable connection to DuckDB
/// * `distances` - Vector of distances (query_name, ref_name, distance) to insert
/// * `level` - Initial level from which to start the search (inclusive)
/// * `ctx` - Context with level and parameter information
/// * `query` - Query object that will be modified with the clique results
///
/// # Process
/// 1. Inserts the distances into the edges table
/// 2. Searches for cliques from the specified level upwards
/// 3. If a clique is found, updates the Query and deletes edges
///
/// # Returns
/// Result with () on success, or Error on failure
pub fn look_for_cliques(
    conn: &mut Connection,
    distances: &[(String, String, f64)],
    level: usize,
    levels: &[(String, f64)],
    click_size: usize,
    query: &mut Query,
) -> Result<(), duckdb::Error> {
    // ...
}

/// Updates the code, code_full and code_state fields of `query` at the `level`
/// based on data from the `code` table.
///
/// - If `level == 0`, gets the maximum `code.L_0` from the `code` table.
/// - If `level > 0`, gets the maximum `code.L_{level}`
///   for rows whose `code_full.L_{level-1}` matches `query.code_full[level-1]`.
/// - Assigns `query.code[level]` and `query.code_full[level]` directly.
/// - Marks `query.code_state[level] = "C"`.
pub fn set_clique_code_for_query(
    conn: &mut Connection,
    query: &mut Query,
    level: usize,
) -> Result<(), duckdb::Error> {
    // ...
}

/// Propagates the values of `query` at the `level` to all `samples` in the clique,
/// updating and inserting into the `code`, `code_full` and `code_state` tables.
///
/// # Arguments
/// * `conn`         – Mutable connection to DuckDB.
/// * `clique_nodes` – Slice of sample IDs to update.
/// * `level`        – Taxonomy level (e.g. 0,1,2…).
/// * `query`        – Query object with values already assigned for that level.
///
/// # Returns
/// Result<(), duckdb::Error>
pub fn update_clique_samples(
    conn: &mut Connection,
    clique_nodes: &[String],
    level: usize,
    query: &Query,
) -> Result<(), duckdb::Error> {
    // ...
}

/// Propagates the values of `query` at the `level` to nodes NOT belonging to the clique,
/// updating and inserting into the `code`, `code_full` and `code_state` tables with state "S".
///
/// # Arguments
/// * `conn`       – Mutable connection to DuckDB.
/// * `all_nodes`  – Slice of all considered sample IDs.
/// * `clique_nodes` – Slice of sample IDs that form the clique.
/// * `level`      – Taxonomy level (e.g. 0,1,2…).
/// * `query`      – Query object with values already assigned for that level.
///
/// # Returns
/// Result<(), duckdb::Error>
pub fn update_non_clique_samples(
    conn: &mut Connection,
    all_nodes: &[String],
    clique_nodes: &[String],
    level: usize,
    query: &Query,
) -> Result<(), duckdb::Error> {
    // ...
}

/// Saves all values from query (code, code_full, code_state)
/// into the code, code_full and code_state tables of DuckDB for the current sample.
///
/// Requires Query to have the fields:
/// - sample_name: String
/// - code: Vec<usize>
/// - code_full: Vec<String>
/// - code_state: Vec<String>
pub fn update_duckdb(
    conn: &mut Connection,
    query: &Query,
) -> DuckResult<()> {
    // ...
}
