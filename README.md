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

**IMPORTANT**: Experimental sources (ATcT through GRI-MECH) with priorities 5-13 are prioritized over theoretical sources (priorities 0-4) when available. The system ensures that the highest quality source is always used for each species.

## Hierarchical Source Selection Details

The hierarchical source selection system follows these key principles:

1. **Experimental Over Theoretical**: Sources with priority 5-13 (experimental) are always preferred over sources with priority 0-4 (theoretical).

2. **Complete Replacement**: Higher priority sources completely replace data from lower priority sources.

3. **Traversal Order**: Sources are traversed in order of decreasing priority, starting with ATCT (13) and ending with basic theoretical methods (0).

4. **Case-Insensitive Matching**: Source names are matched case-insensitively to ensure consistent selection.

5. **Temperature Range Handling**: The system selects the best source that covers the requested temperature range.

6. **Fallback Mechanism**: If no source covers the entire temperature range, the system falls back to the highest priority source available.

The implementation in `src/JThermodynamicsData/utils/json_storage.jl` ensures rigorous adherence to these principles, providing consistent and reliable results.

## Complete Data Pipeline

The package includes a comprehensive pipeline to download and process real-world thermodynamic data from all sources:

### Running the Complete Workflow

```bash
julia run_complete_workflow.jl
```

This single command will:
1. Download data from all original sources (ATcT, Burcat, NIST, JANAF, TDE, ThermoML, etc.)
2. Process all species data following the hierarchical priority
3. Generate plots with uncertainty visualization for all sources
4. Create specialized plots showing only the best source for each species
5. Create the database with all thermodynamic data
6. Verify the hierarchical selection is working correctly
7. Generate comprehensive summary reports

### Individual Pipeline Components

If you need to run only specific parts of the pipeline:

1. **Data Download**: 
   ```bash
   julia scripts/download_all_sources.jl
   ```
   Downloads data from all original sources.

2. **Fetch All Sources**:
   ```bash
   julia scripts/fetch_all_sources.jl
   ```
   Processes data for all species from all sources.

3. **Generate All Plots**:
   ```bash
   julia run_all_species_plots.jl
   ```
   Creates comparison plots for each species showing all sources.

4. **Best Source Plots**:
   ```bash
   julia run_best_source_plots.jl
   ```
   Creates plots showing only the best source for each species with uncertainty.

5. **Hierarchical Verification**:
   ```bash
   julia verify_pipeline.jl
   ```
   Verifies that the hierarchical selection system is working correctly.

6. **Generate Source Summary**:
   ```bash
   julia scripts/generate_source_summary.jl
   ```
   Creates detailed reports on which source was used for each species.

## Real Data Implementation

All data sources now fetch actual experimental data:

1. **Authentic Data**: The system uses data downloaded directly from original sources (ATcT, Burcat, NIST, JANAF, etc.)

2. **Source Metadata**: Each data point includes source attribution, timestamp, and uncertainty estimates

3. **Proper Fallback**: The system falls back to progressively lower priority sources only when necessary

4. **Hierarchical Selection**: Always selects the most accurate source available for each species

5. **Uncertainty Visualization**: Plots show uncertainty bands around property curves

## Key Scripts

- `run_complete_workflow.jl`: Complete data pipeline with all steps
- `run_all_species_plots.jl`: Generate plots showing all sources for each species
- `run_best_source_plots.jl`: Generate plots showing only the best source for each species
- `scripts/download_all_sources.jl`: Download data from all sources
- `scripts/fetch_all_sources.jl`: Process data from all sources for all species
- `scripts/generate_source_summary.jl`: Generate detailed source usage reports
- `verify_pipeline.jl`: Verify the hierarchical selection system

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

# Initialize the package (loads all config files and connects to database)
config, db = JThermodynamicsData.initialize()

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

# Get the best source for a species
best_source = JThermodynamicsData.get_best_source_data("O3")
println("Best source for O3: $(best_source["source"]) (priority: $(best_source["priority"]))")

# Run the full pipeline
include("run_full_pipeline.jl")
```

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