#!/usr/bin/env julia

# Example usage of the JThermodynamicsData package
# This script demonstrates basic queries and plots

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using YAML
using DataFrames
using Plots

println("JThermodynamicsData Example Usage")
println("=================================")

# Load configuration
config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
config = load_config(config_path)

# Initialize database connection
db_path = config["general"]["database_path"]
conn = init_database(config)

# List available species
println("Listing first 10 available species...")
species_list = list_available_species(conn, limit=10)
println(species_list)
println()

# List available databases
println("Listing available thermodynamic databases...")
databases = list_available_databases(conn)
println(databases)
println()

# Example 1: Query a species
println("Example 1: Query information about a species")
println("-------------------------------------------")
species_name = "H2O"
species_info = query_species(conn, species_name, config)
println("Species: $(species_info["name"])")
println("Formula: $(species_info["formula"])")
println("CAS Number: $(species_info["cas_number"])")
println("Molecular Weight: $(species_info["molecular_weight"]) g/mol")
println("Temperature coverage: $(species_info["coverage_percentage"])%")

println("Available data sources:")
for source in species_info["data_sources"]
    println("  - $(source["data_source"]): $(source["temperature_range"][1])-$(source["temperature_range"][2]) K")
end
println()

# Example 2: Query properties at a specific temperature with progressive refinement
println("Example 2: Query properties at a specific temperature with progressive refinement")
println("---------------------------------------------------------------------")
temperature = 1000.0  # K
println("Properties of $species_name at $temperature K:")

# Get all sources to show the refinement process
all_sources = query_properties(conn, species_name, temperature, config, return_all_sources=true)

println("Progressive refinement of data:")
for (i, source) in enumerate(all_sources["results"])
    source_name = source["data_source"]
    println("  Source #$i: $source_name")
    
    for (prop_name, prop_data) in source["properties"]
        println("    $prop_name: $(prop_data["value"]) ± $(prop_data["uncertainty"]) $(prop_data["units"])")
    end
end

# Get final refined data
properties = query_properties(conn, species_name, temperature, config)
println("\nFinal refined properties:")
for (prop_name, prop_data) in properties["properties"]
    println("  $prop_name: $(prop_data["value"]) ± $(prop_data["uncertainty"]) $(prop_data["units"])")
end
println("Data source: $(properties["data_source"])")
println()

# Example 3: Calculate properties over a temperature range
println("Example 3: Calculate properties over a temperature range")
println("-----------------------------------------------------")
temp_range = [300.0, 2000.0]  # K
step = 100.0  # K
println("Calculating properties of $species_name from $(temp_range[1]) to $(temp_range[2]) K:")
range_properties = calculate_properties(conn, species_name, temp_range, config, step=step)

# Plot the data
println("Creating plots...")
p = plot_properties(range_properties, config)
savefig(p, "example_plot.png")
println("Plot saved to 'example_plot.png'")
println()

# Example 4: Compare species
println("Example 4: Compare species")
println("------------------------")
species_list = ["H2O", "CO2", "CH4"]
println("Comparing heat capacity (Cp) for species: $(join(species_list, ", "))")
p_comparison = plot_species_comparison(conn, species_list, "Cp", temp_range, config, step=step)
savefig(p_comparison, "species_comparison.png")
println("Comparison plot saved to 'species_comparison.png'")
println()

# Example 5: Uncertainty analysis
println("Example 5: Uncertainty analysis")
println("-----------------------------")
println("Analyzing uncertainty in enthalpy (H) for $species_name:")
p_uncertainty = plot_uncertainty_analysis(conn, species_name, "H", temp_range, config, step=step)
savefig(p_uncertainty, "uncertainty_analysis.png")
println("Uncertainty analysis saved to 'uncertainty_analysis.png'")
println()

# Example 6: Compare data sources
println("Example 6: Compare data sources")
println("-----------------------------")
println("Comparing data sources for $species_name at $temperature K:")
comparison_table = compare_data_sources_table(conn, species_name, temperature, config)
println(comparison_table)
println()

# Example 7: Theoretical calculation fallback
println("Example 7: Theoretical calculation fallback")
println("----------------------------------------")
theoretical_species = "C6H6"  # Benzene
println("Attempting to get properties for $theoretical_species at 500 K:")
try
    theoretical_properties = query_properties(conn, theoretical_species, 500.0, config)
    println("Source: $(theoretical_properties["data_source"])")
    for (prop_name, prop_data) in theoretical_properties["properties"]
        println("  $prop_name: $(prop_data["value"]) ± $(prop_data["uncertainty"]) $(prop_data["units"])")
    end
catch e
    println("Error: $e")
end
println()

# Close the database connection
close_database(conn)

println("Example script completed successfully!")