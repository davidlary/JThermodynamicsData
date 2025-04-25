# JThermodynamicsData Configuration Guide

This document describes the configuration-driven approach implemented in JThermodynamicsData. All settings that were previously hardcoded are now stored in YAML configuration files, making it easy to customize the behavior of the package without modifying the code.

## Configuration Files

### Main Configuration (`settings.yaml`)

The main configuration file contains general settings for the package, including:

- Database path and settings
- Global temperature range
- Log level
- Data source definitions with priorities
- Theoretical calculation settings

This file serves as the entry point for configuration and defines all data sources with their corresponding priorities.

### Species Configuration (`Species.yaml`)

The species configuration file defines:

- The list of species to process (single source of truth)
- Temperature range for calculations
- Output options for plots and data
- Theoretical calculation options

This file is the single source of truth for which species to process. Ionic species (containing + or -) are automatically detected from this list.

### Advanced Settings (`advanced_settings.yaml`)

This file contains more specific settings that were previously hardcoded:

- **Hierarchical Sources Configuration**:
  - Order of data sources to check (from lowest to highest priority)
  - Fallback settings for when methods fail (default_uncertainty: 1000.0)

- **Uncertainty Configuration**:
  - Monte Carlo samples (1000)
  - Default uncertainties (5% for experimental, 10% for theoretical)
  - How to combine sources (replacement, weighted average)

- **Standard Reference Conditions**:
  - Standard temperature (298.15 K)
  - Standard pressure (101325.0 Pa) 

- **Database Settings**:
  - Cache size
  - Connection timeout
  - Query templates

This file contains settings that previously required code modifications, making it easier to customize the behavior of the package without recompiling.

### Plot Settings (`plots.yaml`)

This file contains all plotting-related settings:

- **General Plot Settings**:
  - Size (1200x900 pixels)
  - DPI (300)
  - Margins (10mm with right margin 20mm for legends)
  - Temperature points (100)

- **Line and Color Settings**:
  - Line width (best source: 2.5, other sources: 1.5)
  - Line style (best source: solid, other sources: dashed)
  - Color strategy (priority-based or distinct colors)

- **Uncertainty Visualization**:
  - Ribbon alpha (0.2)

- **Output Format Options**:
  - File format (png, svg, pdf)
  - Plots directory

This file centralizes all plot styling options that were previously scattered throughout the codebase.

## Using Configuration in Code

JThermodynamicsData provides a simple way to initialize and access these configuration settings:

```julia
using JThermodynamicsData

# Initialize the package
config, db = JThermodynamicsData.initialize()

# Access a configuration value
monte_carlo_samples = config["advanced_settings"]["uncertainty"]["monte_carlo"]["samples"]
println("Using $monte_carlo_samples Monte Carlo samples")

# Access plot settings
plot_size = config["plots_settings"]["general"]["size"]
println("Plots will be size $(plot_size[1])x$(plot_size[2])")

# Standard conditions
std_temp = config["advanced_settings"]["computation"]["reference"]["temperature"]
std_pres = config["advanced_settings"]["computation"]["reference"]["pressure"]
println("Standard conditions: $std_temp K, $std_pres Pa")
```

## Fallback Mechanisms

The code includes fallback mechanisms for missing configuration:

1. If a configuration file is missing, a warning is logged, and the code continues with default values.
2. If an expected setting is missing within a configuration file, a sensible default is used.
3. For critical database operations, multiple fallback paths ensure the code can continue functioning.

## Utility Scripts

Several utility scripts are provided to help with configuration and database management:

- `test_yaml_parsing.jl`: Tests that all YAML files have valid syntax
- `test_config.jl`: Tests loading and using configuration values
- `test_init.jl`: Tests the initialization process
- `db_fix.jl`: Comprehensive script to fix database issues
- `reset_database.jl`: Utility to completely reset the database
- `run_simplified.jl`: Simple workflow to test core functionality

## Troubleshooting

If you encounter issues with configuration:

1. Run `julia test_yaml_parsing.jl` to verify all YAML files have valid syntax
2. Run `julia test_config.jl` to check configuration loading and values
3. If database issues occur, run `julia db_fix.jl` to recreate and fix the database
4. For a complete reset, use `julia reset_database.jl` followed by `julia db_fix.jl`

## Best Practices

1. Always use the `initialize()` function to load configuration, as it handles loading all necessary files
2. Access configuration values through the returned config dictionary
3. When adding new features, use configuration for any values that might need customization
4. Always provide sensible defaults for configuration values in case they're missing
5. Check if a key exists in the configuration before accessing it, or use the `get()` function with a default value

## Reference

See the included configuration files for a complete reference of all available settings:
- `config/settings.yaml`
- `config/Species.yaml`
- `config/advanced_settings.yaml`
- `config/plots.yaml`