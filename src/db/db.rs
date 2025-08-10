/// Structure that encapsulates the DuckDB connection.

    /// Initializes a new DuckDB database at the specified file.

    /// Initializes the `sketches` table.

    /// Adds a Sketch to the sketches table

    /// Reconstructs a SketchManager from the database

    /// Creates the `edges` table with the specified columns if it does not exist.

    // Debug table to store node pairs and distances
    /// Creates the `debug` table with columns `Source`, `Target`, and `dist`.
    /// This table is used to store node pairs and their distances.
    /// The `Source` column is of type VARCHAR, `Target` is VARCHAR,
    /// and `dist` is DOUBLE.

    /// Creates and initializes the `levels` table with data based on the provided Vec<u64>.
    /// The Level values will be "L_0", "L_1", ..., "L_N" where N is the number of elements.
    /// The Dist values correspond to the values of Vec<f64> converted to DOUBLE.

    /// Dynamically creates the `code` table based on the existing levels in the `levels` table.
    /// The table will have a `sample` column (VARCHAR, Primary Key) and columns L_0 to L_N (INTEGER).

    /// Dynamically creates the `code_full` table based on the existing levels in the `levels` table.
    /// The table will have a `sample` column (VARCHAR, Primary Key) and columns L_0 to L_N (VARCHAR).

    /// Dynamically creates the `code_state` table based on the existing levels in the `levels` table.
    /// The table will have a `sample` column (VARCHAR, Primary Key) and columns L_0 to L_N (VARCHAR).

    /// Creates the `metadata` table with information about the experiment.
    /// The table contains metadata of the analysis including genus, acronym, levels, kmer_size, and sketch_size.

    /// Main function that initializes the entire database from a TOML file.
    ///
    /// # Arguments
    ///
    /// * `toml_path` - Path to the TOML file containing the complete configuration
    ///
    /// # Example TOML file
    ///
    /// ```
    /// genus = "Escherichia"
    /// acronym = "EC"
    /// levels = "[0.95,0.98,0.99,0.999,0.9999]"
    /// kmer_size = 21
    /// sketch_size = 1000
    /// ```
    ///
    /// # Functionality
    ///
    /// This function performs the following actions:
    /// 1. Reads and parses the TOML file
    /// 2. Extracts and converts the `levels` field from string to Vec<f64>
    /// 3. Converts the f64 values to u64 by multiplying by 100000 (to preserve precision)
    /// 4. Initializes the `metadata` table with the TOML data
    /// 5. Initializes the `levels` table with the parsed values
    /// 6. Creates all other tables (edges, code, code_full, code_state)

/// Helper function to parse the levels string to Vec<f64>
/// 
/// # Arguments
/// 
/// * `levels_str` - String in format "[0.95,0.97,0.98,0.99,0.999]"
/// 
/// # Returns
/// 
/// anyhow::Result<Vec<f64>> - Vector of f64 values or error

    // Prints the query and parameters for debugging

    /// Reads a TOML file and inserts the data into the metadata table.
    /// This function is independent and can be used separately.
    ///
    /// # Arguments
    ///
    /// * `toml_path` - Path to the TOML file containing the metadata
    ///
    /// # Example TOML file
    ///
    /// ```
    /// genus = "Prueba"
    /// acronym = "PRB"
    /// levels = "[0.95,0.97,0.98,0.99,0.999]"
    /// kmer_size = 21
    /// sketch_size = 1200
    /// click_size = 8
    /// click_threshold = 0.8
    /// reference_size = 100
    /// ```

    /// Access to the internal connection, if required for advanced operations.

/// Inserts a row in the code-like table (code, code_full or code_state)
/// with the sample field equal to `sample_id` and the rest of the columns as NULL.
/// Assumes the name of the ID column is `sample` and the rest of the columns may vary and are after.

/// Copies L_0 to L_L fields from ref_id to input_id in the specified table.

/// Inserts a serialized Sketch object into the sketches table

/// Gets all Sketch objects from the table (to reconstruct SketchManager)

// At the end of your duckdb.rs file

    // Your existing tests...
    // New tests for sketches

        // Verify that the sketches table was created correctly

        // ... test code

        // ... test code
