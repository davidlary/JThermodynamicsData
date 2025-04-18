# JThermodynamicsData

A comprehensive Julia package for handling, processing, and integrating thermodynamic data from multiple sources.

## Overview

JThermodynamicsData provides a unified interface to thermodynamic data for chemical species, combining data from multiple sources with a hierarchical prioritization system. The package includes:

- Parsers for common thermodynamic data formats (NASA-7/9, CHEMKIN, JANAF, ThermoML, Burcat, ATcT, TDE)
- Database integration with DuckDB for efficient storage and retrieval
- Theoretical calculation methods for missing data (group contribution, statistical thermodynamics, etc.)
- Uncertainty quantification and propagation
- Visualization tools for comparing data sources

## Project Structure

- `src/` - Source code for the package
  - `JThermodynamicsData.jl` - Main module definition
  - `core/` - Core functionality (config, constants, types)
  - `database/` - Database interface (connection, schema, queries)
  - `parsers/` - Data format parsers for different sources
  - `models/` - Theoretical calculation models
  - `utils/` - Utility functions (caching, conversion, etc.)
  - `visualization/` - Plotting and comparison tools
  - `api/` - Query interfaces

- `config/` - Configuration files
  - `settings.yaml` - General settings and data source priorities
  - `species.yaml` - Species list and processing options

- `data/` - Data storage
  - `cache/` - Cached data from web sources
  - `external/` - External data files (JANAF, NASA, etc.)
  - `test_data/` - Sample data for testing

- `scripts/` - Utility scripts
  - `initialize_database.jl` - Populate the database
  - `update_database.jl` - Update with new data
  - `example_usage.jl` - Usage examples
  - `debug_scripts/` - Debugging tools

## Setup and Installation

1. Clone the repository:
   ```
   git clone https://github.com/your-username/JThermodynamicsData.git
   cd JThermodynamicsData
   ```

2. Activate the project and install dependencies:
   ```julia
   julia> ]
   pkg> activate .
   pkg> instantiate
   ```

3. Initialize the database:
   ```julia
   include("scripts/initialize_database.jl")
   ```

## Usage

### Basic Usage

```julia
using JThermodynamicsData

# Initialize configuration and database connection
config, conn = initialize()

# Query properties for a species at a specific temperature
properties = query_properties(conn, "N2", 298.15, config)
println("Heat capacity of N2: $(properties["properties"]["Cp"]["value"]) J/mol/K")

# Clean up
close_database(conn)
```

### Processing Multiple Species

```julia
using JThermodynamicsData

# Process thermodynamic data for all species in the configuration
include("run_thermodynamics.jl")
```

### Generating Comprehensive Documentation

The package can automatically generate detailed documentation for each species showing all data sources used in the calculation:

```julia
using JThermodynamicsData

# Initialize configuration and database connection
config, conn = initialize()

# Get a species with all its data sources
result, all_sources = progressively_refine_thermodynamic_data(conn, "H2O", 298.15, config)

# Generate Markdown documentation
markdown_file = create_markdown_documentation("H2O", result, all_sources, "docs")
println("Markdown documentation created at: $markdown_file")

# Generate JSON documentation
json_file = create_json_documentation("H2O", result, all_sources, "docs")
println("JSON documentation created at: $json_file")

# Close connection
close_database(conn)
```

The documentation includes:
- Summary of all sources used (theoretical and experimental)
- Final calculated properties with uncertainties
- Detailed table of all data sources with their property values
- List of missing sources from the hierarchy
- Temperature used for the calculation

### Working with Direct Calculator

For debugging purposes, you can use the direct calculator for N2:

```julia
include("scripts/direct_n2_calculator.jl")
result = calculate_n2_properties(298.15)
println("Heat capacity of N2: $(result["properties"]["Cp"]["value"]) J/mol/K")
```

## Hierarchical Data Processing

JThermodynamicsData uses a hierarchical approach to thermodynamic data with thorough traversal of all available data sources for each species:

1. **Theoretical Methods** (lowest priority, used when no experimental data is available)
   - Enhanced Benson Group Additivity Method
   - Statistical Thermodynamics with Quantum Corrections
   - Joback Method for organic molecules
   - Machine Learning models
   - Quantum Chemistry calculations with anharmonicity corrections

2. **Experimental Databases** (higher priority, in ascending order of accuracy)
   - GRI-MECH (1) - Small set of well-validated combustion species
   - CHEMKIN (2) - Combustion mechanisms
   - NASA-CEA (3) - NASA Chemical Equilibrium with Applications
   - JANAF (4) - NIST-JANAF Thermochemical Tables
   - ThermoML (5) - IUPAC standard for thermodynamic data
   - ThermoData-Engine (6) - NIST TDE with critically evaluated data
   - Burcat (7) - Comprehensive thermochemical database
   - ATcT (8, highest priority) - Active Thermochemical Tables

The system ensures that:
1. Every available source is thoroughly checked for each species
2. Data from higher priority sources refines values from lower priority sources
3. Uncertainty is properly propagated through the hierarchical refinement
4. All data sources used are documented in Markdown and JSON format
5. Log-scale temperature plots provide better visualization across wide temperature ranges

## Debugging

If you encounter issues, check the `DEBUG_REPORT.md` file for common problems and their solutions. Several debugging scripts are provided:

- `scripts/debug_single_species.jl` - Test processing of a single species
- `scripts/direct_n2_calculator.jl` - Direct calculation of N2 properties
- `scripts/robust_debug.jl` - Comprehensive debugging with all fixes applied

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- The thermodynamic data used in this package comes from various public databases maintained by NIST, NASA, and other organizations.
- Thanks to the Julia community for providing excellent libraries for scientific computing.