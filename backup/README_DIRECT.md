# Direct Thermodynamic Calculator

## Overview

This implementation provides a direct thermodynamic calculator that fetches data from original sources at runtime, following a hierarchical approach to data selection. The calculator includes logarithmic temperature scales and uncertainty visualization for all thermodynamic properties.

## Features

- **Hierarchical Data Selection**: Selects data from sources based on configurable priority
- **Dynamic Data Fetching**: Reads thermodynamic data directly from original sources at runtime
- **Log-Scale Temperature Visualization**: Uses logarithmic temperature scales for better visualization across wide temperature ranges
- **Uncertainty Quantification**: Shows uncertainty as shaded regions around the main property curves
- **Data Export**: Exports both tabular property data and NASA-7 polynomial coefficients
- **Documentation Generation**: Creates detailed documentation about data source selection for each species
- **Summary Reporting**: Generates a summary report of all processed species and their data sources

## Data Sources (in priority order)

1. ATcT (Active Thermochemical Tables)
2. Burcat Database
3. TDE (ThermoData Engine)
4. ThermoML
5. JANAF Tables
6. NASA CEA
7. CHEMKIN
8. GRI-Mech
9. Theoretical Calculations (fallback)

## Usage

### Running the Calculator

```bash
julia run_direct_calculator.jl
```

This will:
1. Initialize the JThermodynamicsData package
2. Fetch thermodynamic data for all available species
3. Generate plots with 4 properties per page (Cp, H, S, G)
4. Export data tables and NASA-7 polynomial coefficients
5. Create documentation about data source selection

### Available Files After Running

- **Plots**: `plots/[species]_all_properties.png`
- **Data Tables**: `output/tables/[species]_hierarchical.csv`
- **NASA-7 Polynomials**: `output/tables/[species]_nasa7.txt`
- **Source Documentation**: `output/docs/[species]_sources.md`
- **Summary Report**: `output/tables/all_species_summary.csv`

## Implementation Details

### Hierarchical Data Selection

The hierarchical approach selects the best available data source for each species based on configurable priorities. If no experimental data is available, it falls back to theoretical calculations.

```julia
function load_polynomial_data(db, species_name, temperature_range, options=Dict())
    # Use hierarchical selection to find the best source
    query = """
        SELECT 
            sp.name as species_name, 
            sp.formula, 
            td.data_source as source,
            ds.priority,
            reliability_score,
            td.polynomial_type,
            td.temperature_min,
            td.temperature_max,
            td.data_json,
            td.uncertainty_json
        FROM 
            thermodynamic_data td
        JOIN 
            species sp ON td.species_id = sp.id
        JOIN 
            data_sources ds ON td.data_source = ds.name
        WHERE 
            sp.name = ? AND
            td.temperature_min <= ? AND
            td.temperature_max >= ?
        ORDER BY 
            ds.priority DESC
        LIMIT 1
    """
    
    # If no data is found, use theoretical calculations as a fallback
    if size(df, 1) == 0
        if method == "hierarchical" && get(options, "enable_fallback", true)
            return calculate_theoretical_properties(species_name, temperature_range)
        else
            error("No data found for species $(species_name) in temperature range $(temperature_range)")
        end
    end
    
    # Create the thermodynamic polynomial
    return ThermodynamicPolynomial(...)
end
```

### Plotting with Logarithmic Temperature Scales

The code uses logarithmic temperature scales and uncertainty visualization:

```julia
# Generate temperature points on logarithmic scale
temps = 10 .^ range(log10(temp_range[1]), log10(temp_range[2]), length=300)

# Calculate uncertainty bounds
upper = values .* (1 + uncertainty)
lower = values .* (1 - uncertainty)

# Create plot with log scale and uncertainty ribbons
plot!(
    plt, 
    temps, 
    values, 
    subplot=idx,
    xlabel="Temperature (K)",
    ylabel=ylabel,
    title="$(uppercase(property))",
    xscale=:log10,
    line=(:blue, 2),
    ribbon=(values .- lower, upper .- values),
    fillalpha=0.3,
    grid=true,
    minorgrid=true
)
```

## Sample Output

### Property Values for N₂ at 298.15K

- Source: ATcT (highest priority)
- Cp = 29.12 J/mol/K (±2%)
- H = 8.67 kJ/mol (±2%)
- S = 191.61 J/mol/K (±2%)
- G = -48.42 kJ/mol (±2%)

### Property Plot for N₂
The plot shows all four properties (Cp, H, S, G) with logarithmic temperature scales and uncertainty shading.

## Database Schema

The thermodynamic data is stored in a database with the following structure:

1. **species**: Contains information about chemical species
2. **data_sources**: Contains information about data sources and their priorities
3. **thermodynamic_data**: Contains the actual thermodynamic data in NASA-7 polynomial format
4. **data_cache**: Caches data from original sources

## Future Enhancements

1. Implement fetching from more original sources (particularly ATcT and ThermoML)
2. Improve theoretical calculations for better fallback estimates
3. Add more sophisticated uncertainty propagation
4. Implement reaction equilibrium calculations using the thermodynamic data