#!/usr/bin/env julia

# A simplified script that uses the direct N2 calculator to process N2
# This script bypasses all the complex database machinery and uses our standalone implementation

using Pkg
Pkg.activate(dirname(@__DIR__))

using JThermodynamicsData
using YAML
using Plots
using DataFrames

println("JThermodynamicsData - Direct N2 Processor")
println("======================================")

# Load species configuration 
species_config_path = joinpath(dirname(@__DIR__), "config", "species.yaml")
species_config = YAML.load_file(species_config_path)

# Extract settings from species config
temp_range = get(species_config, "temperature_range", [100.0, 10000.0])
temp_step = get(species_config, "temperature_step", 10.0)
output_options = get(species_config, "output", Dict())

# Create output directories
plot_dir = get(output_options, "plot_directory", "plots")
mkpath(plot_dir)
tables_dir = joinpath(dirname(@__DIR__), "output", "tables")
mkpath(tables_dir)

# Include the direct calculator
include(joinpath(dirname(@__DIR__), "scripts", "direct_n2_calculator.jl"))

println("Processing N2 with direct calculator")
println("Temperature range: $(temp_range[1])-$(temp_range[2]) K")
println("Step size: $temp_step K")

# Calculate properties across the temperature range
result = calculate_n2_properties_range(temp_range, temp_step)

# Generate plots
println("Generating plots...")

# Heat capacity plot
temps = result["temperatures"]
cp_values = result["properties"]["Cp"]["values"]
cp_uncertainties = result["properties"]["Cp"]["uncertainties"]

p_cp = plot(
    temps,
    cp_values,
    ribbon=cp_uncertainties,
    fillalpha=0.3,
    xlabel="Temperature (K)",
    ylabel="Cp (J/mol/K)",
    title="Heat Capacity for N2",
    legend=false,
    lw=2,
    grid=true
)

# Enthalpy plot
h_values = result["properties"]["H"]["values"]
h_uncertainties = result["properties"]["H"]["uncertainties"]

p_h = plot(
    temps,
    h_values,
    ribbon=h_uncertainties,
    fillalpha=0.3,
    xlabel="Temperature (K)",
    ylabel="H (kJ/mol)",
    title="Enthalpy for N2",
    legend=false,
    lw=2,
    grid=true
)

# Entropy plot
s_values = result["properties"]["S"]["values"]
s_uncertainties = result["properties"]["S"]["uncertainties"]

p_s = plot(
    temps,
    s_values,
    ribbon=s_uncertainties,
    fillalpha=0.3,
    xlabel="Temperature (K)",
    ylabel="S (J/mol/K)",
    title="Entropy for N2",
    legend=false,
    lw=2,
    grid=true
)

# Gibbs energy plot
g_values = result["properties"]["G"]["values"]
g_uncertainties = result["properties"]["G"]["uncertainties"]

p_g = plot(
    temps,
    g_values,
    ribbon=g_uncertainties,
    fillalpha=0.3,
    xlabel="Temperature (K)",
    ylabel="G (kJ/mol)",
    title="Gibbs Energy for N2",
    legend=false,
    lw=2,
    grid=true
)

# Combine into a 2x2 layout
p = plot(p_cp, p_h, p_s, p_g, layout=(2, 2), size=(1000, 800))

# Save the plot
plot_format = get(output_options, "plot_format", "png")
plot_file = joinpath(plot_dir, "N2_properties_direct.$plot_format")
savefig(p, plot_file)
println("Saved properties plot to: $plot_file")

# Create tabular output
println("Creating tabular output...")
table_file = joinpath(tables_dir, "N2_data_direct.csv")

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
println("Direct processing complete!")