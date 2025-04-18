# JThermodynamicsData

A comprehensive Julia package for handling, processing, and integrating thermodynamic data from multiple sources.

## Overview

JThermodynamicsData provides a unified interface to thermodynamic data for chemical species, combining data from multiple sources with a hierarchical prioritization system. The package includes:

- Parsers for common thermodynamic data formats (NASA-7/9, CHEMKIN, JANAF, ThermoML, Burcat, ATcT, TDE)
- Database integration with DuckDB for efficient storage and retrieval
- Theoretical calculation methods for missing data (group contribution, statistical thermodynamics, etc.)
- Uncertainty quantification and propagation
- Visualization tools for comparing data sources
- CAS Registry Number integration for reliable species identification

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
  - `test_data/` - Sample data for testing

- `scripts/` - Utility scripts
  - `initialize_database.jl` - Populate the database
  - `update_database.jl` - Update with new data
  - `example_usage.jl` - Usage examples

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

## Hierarchical Thermodynamic Data Sources

JThermodynamicsData implements a comprehensive hierarchical approach to thermodynamic data processing. The system progressively refines property estimates by incorporating data from sources in order of increasing reliability. All theoretical approaches are used and displayed in visualization, but for species with experimental data, only experimental values are used in the final calculation.

### Table of Data Sources

| Priority | Source Type             | Source                  | Reliability (1-5) | Description                                          | Citation                                                                  |
|----------|-------------------------|-------------------------|-------------------|------------------------------------------------------|---------------------------------------------------------------------------|
| 0        | Theoretical (Level 1)   | Group Contribution      | 1.0               | Basic group contribution methods                     | Joback & Reid (1987) [DOI: 10.1002/ceat.270100209]                        |
| 0        | Theoretical (Level 2)   | Statistical Thermodynamics | 2.0           | Basic statistical mechanics                          | McQuarrie (2000) Statistical Mechanics                                    |
| 0        | Theoretical (Level 3)   | Benson Group Additivity | 3.0              | Advanced group contribution                          | Benson et al. (1969) [DOI: 10.1021/cr60259a002]                           |
| 0        | Theoretical (Level 4)   | Quantum-Statistical     | 4.0              | Quantum-corrected statistical mechanics              | Barone (2004) [DOI: 10.1021/jp049822w]                                    |
| 1        | Experimental (Level 1)  | GRI-MECH                | 3.0              | Gas Research Institute mechanism                     | Smith et al. (1999) [DOI: 10.1.1.16.309]                                  |
| 2        | Experimental (Level 2)  | CHEMKIN                 | 3.5              | Combustion chemistry mechanism                       | Kee et al. (1996) CHEMKIN-III                                             |
| 3        | Experimental (Level 3)  | NASA-CEA                | 4.0              | NASA Chemical Equilibrium with Applications          | Gordon & McBride (1994) NASA RP-1311                                      |
| 4        | Experimental (Level 4)  | JANAF                   | 4.5              | NIST-JANAF Thermochemical Tables                     | Chase (1998) [DOI: 10.18434/T42S31]                                       |
| 5        | Experimental (Level 5)  | ThermoML                | 4.5              | IUPAC standard for thermodynamic data                | Frenkel et al. (2011) [DOI: 10.1063/1.3525836]                            |
| 6        | Experimental (Level 6)  | ThermoData Engine (TDE) | 4.8              | NIST TDE with critically evaluated data              | Diky et al. (2012) [DOI: 10.1021/je300128w]                               |
| 7        | Experimental (Level 7)  | Burcat Database         | 4.9              | Burcat & Ruscic thermochemical database              | Burcat & Ruscic (2005) [DOI: 10.2172/925269]                              |
| 8        | Experimental (Level 8)  | Active Thermo Tables (ATcT) | 5.0         | Active Thermochemical Tables with network approach   | Ruscic et al. (2005) [DOI: 10.1063/1.1804602]                             |

### Data Source Hierarchy and Reliability

The hierarchical approach ensures that:

1. **Progressive Refinement**: Data is progressively refined from theoretical approaches to increasingly accurate experimental sources.
2. **Theory vs. Experiment**: All theoretical approaches are displayed in visualization, but for species with experimental data, only experimental values contribute to the final calculation.
3. **Uncertainty Propagation**: Uncertainty is properly propagated through the hierarchical refinement process.
4. **Documentation**: All data sources used for each species are thoroughly documented.

#### Theoretical Methods (In Order of Increasing Accuracy)

1. **Group Contribution Methods** (Reliability: Low)
   - Simple additive models based on molecular fragments
   - Uses Joback & Reid method for small organic molecules
   - Fastest but least accurate approach
   - Typical uncertainty: 10-30%

2. **Statistical Thermodynamics** (Reliability: Medium-Low)
   - Basic statistical mechanical approach with rigid-rotor and harmonic oscillator approximations
   - Requires molecular parameters (geometries, vibrational frequencies)
   - More accurate than group contribution but with significant approximations
   - Typical uncertainty: 5-15%

3. **Benson Group Additivity** (Reliability: Medium)
   - Advanced group contribution method with next-nearest neighbor corrections
   - Accounts for ring corrections and other non-additive effects
   - More accurate than basic statistical thermodynamics for many species
   - Typical uncertainty: 5-10%

4. **Quantum-Statistical Thermodynamics** (Reliability: Medium-High)
   - Quantum-corrected statistical mechanics with anharmonicity effects
   - Includes internal rotation treatments, anharmonic corrections
   - Most accurate theoretical approach implemented in the package
   - Typical uncertainty: 3-8%

#### Experimental Databases (In Order of Increasing Reliability)

Experimental sources are assigned priority levels from 1 (lowest) to 8 (highest), determining the order in which they are consulted:

1. **GRI-MECH 3.0** (Priority: 1)
   - Focused on combustion chemistry for natural gas
   - Limited set of species (~50) but well-validated
   - Used primarily for combustion modeling

2. **CHEMKIN Mechanisms** (Priority: 2)
   - Various reaction mechanisms in CHEMKIN format
   - Coverage depends on the specific mechanism
   - Often based on older NASA polynomial fits

3. **NASA Chemical Equilibrium with Applications (CEA)** (Priority: 3)
   - NASA's thermodynamic database
   - Covers ~2000 species
   - Good for high-temperature applications

4. **NIST-JANAF Thermochemical Tables** (Priority: 4)
   - Gold standard for many years
   - Well-documented uncertainty estimates
   - Updated periodically by NIST

5. **ThermoML** (Priority: 5)
   - IUPAC standard for thermodynamic data exchange
   - Includes experimental metadata
   - Covers many complex organic compounds

6. **NIST ThermoData Engine (TDE)** (Priority: 6)
   - Critically evaluated thermodynamic data
   - Includes uncertainty quantification
   - Covers thousands of organic compounds

7. **Burcat Database** (Priority: 7)
   - Comprehensive database maintained by Burcat & Ruscic
   - Regular updates with new species
   - Focus on species relevant to combustion and atmospheric chemistry

8. **Active Thermochemical Tables (ATcT)** (Priority: 8)
   - Most accurate source, using thermochemical network approach
   - Self-consistent thermodynamics across all species
   - Rigorous uncertainty quantification
   - Limited species coverage but growing

### Output Formats

JThermodynamicsData generates comprehensive output for each species:

1. **NASA-7 and NASA-9 Polynomials**
   - Coefficients for both formats
   - Uncertainty polynomials
   - Temperature-specific uncertainty estimates

2. **CAS Registry Numbers**
   - Included in all outputs for reliable species identification
   - Cross-referenced across all data sources

3. **Comprehensive Documentation**
   - Markdown and JSON formats
   - Source attribution
   - Uncertainty documentation
   - Calculation methodology

4. **Data Tables**
   - CSV spreadsheets with data across temperature range
   - Property values and uncertainties
   - Units explicitly specified

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

## Working with Direct Calculator

For debugging purposes, you can use the direct calculator:

```julia
include("hierarchical_calculator.jl")
result = calculate_properties("N2", 298.15)
println("Heat capacity of N2: $(result["properties"]["Cp"]["value"]) J/mol/K")
```

## Data Processing Workflow

The data processing workflow follows these steps:

1. **Data Source Retrieval**
   - Retrieve data from original sources
   - Implement local caching with version checking
   - Check for newer versions when cached files exist

2. **Theoretical Calculations**
   - Apply theoretical methods in order of increasing accuracy
   - Store uncertainties for each method
   - Display all methods in visualizations

3. **Experimental Data Processing**
   - Process all available experimental data sources in order of priority
   - Refine results progressively with higher-quality sources
   - Propagate uncertainties through the refinement process

4. **Output Generation**
   - Generate NASA-7 and NASA-9 polynomial fits
   - Create uncertainty polynomials
   - Produce documentation in multiple formats
   - Generate CSV data tables with uncertainties

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