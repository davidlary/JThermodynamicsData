# JThermodynamicsData Scripts

This directory contains utility scripts for working with the JThermodynamicsData package.

## Main Utility Scripts

- `initialize_database.jl` - Initialize the database with data from all configured sources
- `update_database.jl` - Update the database with new data
- `example_usage.jl` - Examples of using the JThermodynamicsData API

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

## Notes on Fixes

The debugging scripts implement several important fixes that may need to be incorporated into the main codebase:

1. Constant name matching (lowercase vs uppercase)
2. Type conversion handling (Real to Float64)
3. Proper error handling for missing data
4. Explicit keyword argument passing
5. Safe metadata parsing

For full details on all issues and fixes, see the `DEBUG_REPORT.md` file in the main directory.