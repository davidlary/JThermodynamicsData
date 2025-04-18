#!/usr/bin/env julia

# Debug script to process just the first species (N2) from the species list
# This script is for debugging purposes to verify the pipeline works without errors

using Pkg
Pkg.activate(".")

println("Loading JThermodynamicsData module...")
using JThermodynamicsData
using YAML
using Plots
using Dates

# Define the log-temperature plotting function
function plot_properties_log_temp(property_data::Dict, config::Dict)
    # Create individual plots
    temps = property_data["temperatures"]
    
    # Extract property data
    cp_values = property_data["properties"]["Cp"]["values"]
    cp_uncertainties = property_data["properties"]["Cp"]["uncertainties"]
    
    h_values = property_data["properties"]["H"]["values"]
    h_uncertainties = property_data["properties"]["H"]["uncertainties"]
    
    s_values = property_data["properties"]["S"]["values"]
    s_uncertainties = property_data["properties"]["S"]["uncertainties"]
    
    g_values = property_data["properties"]["G"]["values"]
    g_uncertainties = property_data["properties"]["G"]["uncertainties"]
    
    # Set up plots
    species_name = property_data["species_name"]
    
    # Heat capacity plot
    p_cp = plot(
        temps,
        cp_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="Cp (J/mol/K)",
        title="Heat Capacity for $species_name",
        ribbon=cp_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Enthalpy plot
    p_h = plot(
        temps,
        h_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="H (kJ/mol)",
        title="Enthalpy for $species_name",
        ribbon=h_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Entropy plot
    p_s = plot(
        temps,
        s_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="S (J/mol/K)",
        title="Entropy for $species_name",
        ribbon=s_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Gibbs energy plot
    p_g = plot(
        temps,
        g_values,
        xscale=:log10,
        xlabel="Temperature (K)",
        ylabel="G (kJ/mol)",
        title="Gibbs Energy for $species_name",
        ribbon=g_uncertainties,
        fillalpha=0.3,
        legend=false,
        lw=2,
        grid=true,
        framestyle=:box
    )
    
    # Combine into a 2x2 layout
    p = plot(p_cp, p_h, p_s, p_g, layout=(2, 2), size=(1000, 800))
    
    return p
end

println("Start time: $(now())")

# Load config and species
println("Loading configuration...")
config_path = joinpath(dirname(@__DIR__), "config", "settings.yaml")
config = JThermodynamicsData.load_config(config_path)

# Initialize database connection
println("Initializing database...")
conn = JThermodynamicsData.init_database(config)

# Load species list but override to only include the first species (N2)
species_config_path = joinpath(dirname(@__DIR__), "config", "species.yaml")
species_config = YAML.load_file(species_config_path)

# Extract just the first species
species_name = species_config["species"][1]
println("Processing single species for debugging: $species_name")

# Extract settings from species config
temp_range = get(species_config, "temperature_range", [100.0, 10000.0])
temp_step = get(species_config, "temperature_step", 10.0)
output_options = get(species_config, "output", Dict())

# Create output directories
if get(output_options, "generate_plots", false)
    plot_dir = get(output_options, "plot_directory", "plots")
    mkpath(plot_dir)
end

tables_dir = joinpath(dirname(@__DIR__), "output", "tables")
mkpath(tables_dir)

try
    # Check if species exists in the database
    println("Checking if species exists in database...")
    species_id = JThermodynamicsData.get_species_id(conn, species_name)
    
    if species_id === nothing
        println("Species not found in database. Will use theoretical calculations only.")
    else
        println("Species found in database with ID: $species_id")
    end
    
    # Calculate properties over the temperature range
    println("Calculating thermodynamic properties...")
    result = JThermodynamicsData.calculate_properties(conn, species_name, temp_range, config; step=temp_step)
    
    println("Temperature points calculated: $(length(result["temperatures"]))")
    
    # Generate plots if requested
    if get(output_options, "generate_plots", true)
        use_log_temp = get(output_options, "use_log_temperature", false)
        plot_format = get(output_options, "plot_format", "png")
        plot_dir = get(output_options, "plot_directory", "plots")
        
        println("Generating plots...")
        
        # Create plots with appropriate temperature scale
        if use_log_temp
            p = plot_properties_log_temp(result, config)
        else
            p = JThermodynamicsData.plot_properties(result, config)
        end
        
        # Save the plot
        plot_file = joinpath(plot_dir, "$(species_name)_properties.$plot_format")
        savefig(p, plot_file)
        println("Saved properties plot to: $plot_file")
        
        # Create uncertainty analysis plot
        p_uncertainty = JThermodynamicsData.plot_uncertainty_analysis(conn, species_name, "Cp", temp_range, config, step=temp_step)
        uncertainty_file = joinpath(plot_dir, "$(species_name)_uncertainty.$plot_format")
        savefig(p_uncertainty, uncertainty_file)
        println("Saved uncertainty plot to: $uncertainty_file")
        
        # Create data source comparison if available
        try
            p_comparison = JThermodynamicsData.plot_data_source_comparison(conn, species_name, "Cp", temp_range, config, step=temp_step)
            comparison_file = joinpath(plot_dir, "$(species_name)_sources.$plot_format")
            savefig(p_comparison, comparison_file)
            println("Saved data source comparison plot to: $comparison_file")
        catch e
            println("Could not create data source comparison: $e")
        end
    end
    
    # Create tabular output
    println("Creating tabular output...")
    table_file = joinpath(tables_dir, "$(species_name)_data.csv")
    
    # Open CSV file
    open(table_file, "w") do io
        # Write header
        println(io, "Temperature,Cp,Cp_uncertainty,H,H_uncertainty,S,S_uncertainty,G,G_uncertainty")
        
        # Write data rows
        for i in 1:length(result["temperatures"])
            temp = result["temperatures"][i]
            cp = result["properties"]["Cp"]["values"][i]
            cp_unc = result["properties"]["Cp"]["uncertainties"][i]
            h = result["properties"]["H"]["values"][i]
            h_unc = result["properties"]["H"]["uncertainties"][i]
            s = result["properties"]["S"]["values"][i]
            s_unc = result["properties"]["S"]["uncertainties"][i]
            g = result["properties"]["G"]["values"][i]
            g_unc = result["properties"]["G"]["uncertainties"][i]
            
            println(io, "$temp,$cp,$cp_unc,$h,$h_unc,$s,$s_unc,$g,$g_unc")
        end
    end
    println("Saved tabular data to: $table_file")
    
catch e
    println("ERROR processing $species_name: $e")
    println(stacktrace(catch_backtrace()))
finally
    # Close database connection
    println("Closing database connection...")
    JThermodynamicsData.close_database(conn)
end

println("End time: $(now())")
println("Debug processing complete!")