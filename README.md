# JThermodynamicsData

A comprehensive Julia package for accessing, processing, and visualizing thermodynamic data from multiple sources with built-in hierarchical selection.

## Overview

JThermodynamicsData provides a unified interface to thermodynamic data for chemical species, combining data from multiple sources with a hierarchical prioritization system. The package implements a rigorously structured hierarchy that automatically selects the most accurate source available for each species based on data quality and reliability.

The system is designed to ensure **complete replacement** of data from lower priority sources with higher priority sources, guaranteeing that the final result exactly matches the highest priority source available for each species.

## Key Features

- **Multi-source data access**: Access thermodynamic data from 14 different sources, including:
  - Active Thermochemical Tables (ATcT)
  - Burcat Database
  - NIST Chemistry WebBook
  - NIST ThermoData Engine (TDE)
  - JANAF Thermochemical Tables
  - NASA Chemical Equilibrium with Applications (CEA)
  - And more...

- **Hierarchical selection**: Automatically selects the highest quality data source available for each species.

- **Data visualization**: Generate comparative plots of thermodynamic properties across multiple sources.

- **Database integration**: Store and query thermodynamic data using DuckDB.

- **Comprehensive validation**: Ensure proper hierarchy respect and data accuracy.

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
  - `Species.yaml` - Species list and processing options

- `data/` - Data storage
  - `cache/` - Cached data from web sources
  - `external/` - External data files (JANAF, NASA, etc.)
  - `species/` - JSON storage for species data
  - `test_data/` - Sample data for testing

- `scripts/` - Utility scripts
  - `download_all_sources.jl` - Download data from all sources
  - `fetch_all_sources.jl` - Process and parse data
  - `sync_json_to_database.jl` - Sync data to DuckDB

## Data Hierarchy

The package implements a strict hierarchy of data sources based on their accuracy and reliability:

| Priority | Source | Reliability | Description |
|---------|--------|-------------|-------------|
| 13 | atct | 5.0 | Active Thermochemical Tables (highest accuracy) |
| 12 | burcat | 4.9 | Burcat Database (very high accuracy) |
| 11 | nist-webbook | 4.95 | NIST Chemistry WebBook (high accuracy) |
| 10 | tde | 4.8 | NIST ThermoData Engine (high accuracy) |
| 9 | thermoml | 4.5 | ThermoML Standard (good accuracy) |
| 8 | janaf | 4.5 | JANAF Thermochemical Tables (good accuracy) |
| 7 | nasa-cea | 4.0 | NASA CEA Database (decent accuracy) |
| 6 | chemkin | 3.8 | CHEMKIN format data (decent accuracy) |
| 5 | gri-mech | 3.5 | GRI-MECH 3.0 (decent accuracy) |
| 4 | quantum-statistical | 4.0 | Quantum-Statistical calculations (best theoretical) |
| 3 | benson-group | 3.0 | Benson Group Additivity (good theoretical) |
| 2 | stat-thermo | 3.5 | Statistical thermodynamics (moderate theoretical) |
| 1 | group-contribution | 3.0 | Group contribution methods (basic theoretical) |
| 0 | theoretical | 2.5 | Basic theoretical calculations (lowest accuracy) |

**IMPORTANT**: Experimental sources (ATcT through GRI-MECH) are prioritized over theoretical sources when available. The system ensures that the highest quality source is always used for each species.

## Setup and Installation

1. Clone the repository:
   ```
   git clone https://github.com/your-username/JThermodynamicsData.git
   cd JThermodynamicsData
   ```

2. Set up the environment and install all required packages:
   ```bash
   julia scripts/setup_packages.jl
   ```

3. Run the complete workflow script:
   ```bash
   julia run_complete_workflow.jl
   ```

The complete workflow script handles:
1. Downloading all data sources
2. Processing data into JSON files
3. Syncing data to DuckDB database
4. Generating comparison plots
5. Validating the hierarchy selection

## Workflow

The new consolidated workflow simplifies the thermodynamic data processing:

```bash
# Run the complete workflow (download, process, database, plots, validation)
julia run_complete_workflow.jl
```

This single command replaces the previous multi-step process:
```bash
# Old workflow (no longer needed)
# julia scripts/download_all_sources.jl
# julia scripts/fetch_all_sources.jl
# julia run_all_species_plots.jl
# julia scripts/sync_json_to_database.jl
```

## Troubleshooting

If you encounter any issues:

1. Clear compiled code:
   ```bash
   rm -rf ~/.julia/compiled/v1.*/JThermodynamicsData
   ```

2. Reinitialize the database:
   ```bash
   rm -f ~/Dropbox/Environments/Code/JThermodynamicsData/data/thermodynamics.duckdb
   ```

3. Run the complete workflow again:
   ```bash
   julia run_complete_workflow.jl
   ```

## Hierarchical Thermodynamic Data Sources

JThermodynamicsData implements a comprehensive hierarchical approach to thermodynamic data processing. The system ensures that for each species, the most accurate available source is used, following these principles:

1. **Full Hierarchical Traversal**: The system always processes all sources in order of increasing priority.
2. **Complete Replacement**: Higher priority sources completely replace data from lower priority sources (no weighted averaging).
3. **Most Accurate Source**: The final refined result is exactly the value from the most accurate source available (highest priority).
4. **Theoretical Source Display**: Plots show only ONE theoretical source (best available) labeled consistently as "Theoretical".
5. **Consistent Source Naming**: All sources in plots are displayed with proper capitalization.
6. **Documentation**: All data sources used for each species are thoroughly documented, including which source was selected as the most accurate.

### Theoretical vs. Experimental Sources

The package provides both theoretical calculations and experimental data:

- **Theoretical Methods** (Priority 0-4)
  - Basic theoretical calculations
  - Group contribution methods
  - Statistical thermodynamics
  - Benson Group Additivity
  - Quantum-Statistical Thermodynamics

- **Experimental Databases** (Priority 5-13)
  - GRI-MECH 3.0
  - CHEMKIN format data
  - NASA CEA Database
  - JANAF Thermochemical Tables
  - ThermoML Standard
  - NIST ThermoData Engine
  - NIST Chemistry WebBook
  - Burcat Database
  - Active Thermochemical Tables

The hierarchical system automatically uses experimental data when available, falling back to theoretical calculations only when necessary. For ALL species (including ions like OH-, NO+, NO2+, Zn), the system ensures proper hierarchical traversal.

## Output

The package produces several outputs:

- **Plots**: Saved in `plots/` directory - Shows all available sources with the best one highlighted
- **Data tables**: Saved in `output/tables/` directory - Contains numerical data for all properties
- **Source documentation**: Saved in `output/docs/` directory - Documents sources used for each species
- **Summary report**: Saved to `output/summary_report.md` - Overview of all processed species
- **Database**: Saved to `data/thermodynamics.duckdb` - For efficient querying

## Basic Usage

```julia
using JThermodynamicsData

# Run the complete workflow
include("run_complete_workflow.jl")

# Calculate thermodynamic properties for a species
properties = JThermodynamicsData.calculate_properties("H2O", 298.15)
println("Cp of H2O at 298.15 K: $(properties.cp) J/mol·K")
println("Source used: $(properties.source)")

# See all available sources for a species
sources = JThermodynamicsData.list_available_sources("O3")
println("Available sources for O3:")
for (source, priority, reliability) in sources
    println("  $source (priority: $priority, reliability: $reliability)")
end
```

## Storage Systems

JThermodynamicsData now uses a combined approach:

1. **JSON Storage**: Each species has a JSON file in `data/species/` containing data from all sources
2. **Database Storage**: Data is synced to DuckDB for efficient querying

This combined approach provides both human-readable storage and efficient query capabilities.

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

## References

### Theoretical Methods

1. Joback, K. G., & Reid, R. C. (1987). Estimation of pure-component properties from group-contributions. Chemical Engineering Communications, 57(1-6), 233-243. DOI: 10.1002/ceat.270100209

2. McQuarrie, D. A. (2000). Statistical Mechanics. University Science Books.

3. Benson, S. W., Cruickshank, F. R., Golden, D. M., Haugen, G. R., O'neal, H. E., Rodgers, A. S., ... & Walsh, R. (1969). Additivity rules for the estimation of thermochemical properties. Chemical Reviews, 69(3), 279-324. DOI: 10.1021/cr60259a002

4. Barone, V. (2004). Vibrational zero-point energies and thermodynamic functions beyond the harmonic approximation. The Journal of Chemical Physics, 120(7), 3059-3065. DOI: 10.1021/jp049822w

### Experimental Databases

5. Smith, G. P., Golden, D. M., Frenklach, M., Moriarty, N. W., Eiteneer, B., Goldenberg, M., ... & Qin, Z. (1999). GRI-Mech 3.0. Available at: http://www.me.berkeley.edu/gri_mech/

6. Kee, R. J., Rupley, F. M., & Miller, J. A. (1996). CHEMKIN-III: A FORTRAN chemical kinetics package for the analysis of gas-phase chemical and plasma kinetics. Sandia National Laboratories Report SAND96-8216.

7. Gordon, S., & McBride, B. J. (1994). Computer program for calculation of complex chemical equilibrium compositions and applications. NASA Reference Publication 1311.

8. Chase, M. W. (1998). NIST-JANAF Thermochemical Tables, 4th Edition. Journal of Physical and Chemical Reference Data, Monograph 9, 1-1951. DOI: 10.18434/T42S31

9. Frenkel, M., Chirico, R. D., Diky, V., Yan, X., Dong, Q., & Muzny, C. (2011). ThermoData Engine (TDE): Software implementation of the dynamic data evaluation concept. Journal of Chemical Information and Modeling, 51(10), 2538-2548. DOI: 10.1063/1.3525836

10. Diky, V., Chirico, R. D., Kazakov, A. F., Muzny, C. D., & Frenkel, M. (2012). ThermoData Engine (TDE): Software implementation of the dynamic data evaluation concept. 9. Extensible thermodynamic constraints for pure compounds and new model developments. Journal of Chemical & Engineering Data, 57(10), 2736-2746. DOI: 10.1021/je300128w

11. Burcat, A., & Ruscic, B. (2005). Third millennium ideal gas and condensed phase thermochemical database for combustion with updates from active thermochemical tables. Argonne National Laboratory Argonne, IL. DOI: 10.2172/925269

12. Ruscic, B., Pinzon, R. E., Morton, M. L., von Laszevski, G., Bittner, S. J., Nijsure, S. G., ... & Wagner, A. F. (2005). Introduction to active thermochemical tables: Several "key" enthalpies of formation revisited. The Journal of Physical Chemistry A, 109(42), 9396-9409. DOI: 10.1063/1.1804602