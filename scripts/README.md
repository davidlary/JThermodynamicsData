# JThermodynamicsData Scripts

This directory contains utility scripts for working with the JThermodynamicsData package.

## Main Utility Scripts

- `initialize_database.jl` - Initialize the database with data from all configured sources
- `update_database.jl` - Update the database with new data
- `example_usage.jl` - Examples of using the JThermodynamicsData API
- `generate_species_data.jl` - Generate accurate thermodynamic data for species_data.yaml from first principles

## Debugging Scripts

A set of debugging scripts has been created to help troubleshoot issues with the package.

### Basic Scripts

- `debug_single_species.jl` - Process a single species (N2) to verify the main pipeline works
- `simple_debug.jl` - A simplified version that directly accesses the core functions

### Standalone Calculators

- `direct_n2_calculator.jl` - Standalone calculator for N2 using NASA-7 polynomials
  - Can be used independently of the main package
  - Provides a reference implementation for comparison
  - Generates CSV output for verification

- `evaluate_nasa7.jl` - Direct evaluation of NASA-7 polynomial for N2
  - Simple implementation for direct verification
  - Shows values at standard temperature (298.15 K)

- `fix_statistical_thermodynamics.jl` - Standalone statistical thermodynamics calculator
  - Implements the statistical mechanics approach for N2
  - Provides alternative calculation method for comparison

### Comprehensive Debugging

- `robust_debug.jl` - A comprehensive script that:
  - Applies fixes to known issues
  - Performs step-by-step testing
  - Verifies database structure and content
  - Generates verification plots

- `initialize_minimal_database.jl` - Creates a minimal database with just N2
  - Useful for isolated testing
  - Ensures database structure is correct

## Usage

Most scripts can be run directly from Julia:

```
julia scripts/direct_n2_calculator.jl
```

For debugging purposes, it's recommended to follow this sequence:

1. Initialize a minimal database:
   ```
   julia scripts/initialize_minimal_database.jl
   ```

2. Run the direct calculator for reference values:
   ```
   julia scripts/direct_n2_calculator.jl
   ```

3. Run the robust debug script to test all functionality:
   ```
   julia scripts/robust_debug.jl
   ```

## Modifying the Scripts

Feel free to modify these scripts for your own debugging needs. Some useful modifications:

- Change the test species from N2 to another species
- Adjust temperature ranges to test specific regions
- Change data sources to test different hierarchical combinations
- Add custom visualization of results

## generate_species_data.jl Details

This script generates thermodynamic data for ionic species and other special cases from first principles and established databases.

### Purpose

Instead of hardcoding thermodynamic data in the main code, this script:

1. Fetches data from authoritative sources like NIST Chemistry WebBook and JANAF Tables
2. Performs quantum chemistry calculations (or uses published results) for species without experimental data
3. Applies statistical thermodynamics for well-understood cases like monoatomic ions
4. Generates NASA-7 polynomial coefficients from the calculated properties
5. Exports all data to `../config/species_data.yaml` in the required format

### Data Sources

- **NIST Chemistry WebBook**: Primary source for standard enthalpies, entropies, and ionization potentials
- **JANAF Thermochemical Tables**: Classic reference for thermodynamic data
- **Quantum Chemistry Calculations**: G4, CBS-QB3, and W1 level calculations for ions without experimental data
- **Statistical Thermodynamics**: For simple cases where partition functions can be reliably calculated

### Usage

```bash
julia generate_species_data.jl
```

## Notes on Fixes

The debugging scripts implement several important fixes that may need to be incorporated into the main codebase:

1. Constant name matching (lowercase vs uppercase)
2. Type conversion handling (Real to Float64)
3. Proper error handling for missing data
4. Explicit keyword argument passing
5. Safe metadata parsing

For full details on all issues and fixes, see the `DEBUG_REPORT.md` file in the main directory.