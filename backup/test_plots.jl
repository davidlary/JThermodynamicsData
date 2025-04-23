# Simple script to test ionic species plotting

using Plots
using YAML

include("hierarchical_calculator.jl")

# Configure plot output
output_dir = joinpath(@__DIR__, "plots")
mkpath(output_dir)

# Define species to plot
species_list = ["He+", "Xe", "Xe+"]

# Generate and save plots for each species
for species_name in species_list
    println("Generating plots for $species_name...")
    
    # Run the hierarchical calculation
    result = calculate_properties_range(
        species_name, 
        species_name,  # Using species name as formula for simplicity
        [200.0, 6000.0],  # Temperature range
        100.0,  # Temperature step
        Dict()  # Empty data source dict to use defaults
    )
    
    # Plot and save all properties
    properties = ["Cp", "H", "S", "G"]
    for prop in properties
        println("  Plotting $prop...")
        p = plot_thermodynamic_properties(result, prop)
        
        # Save the plot
        output_file = joinpath(output_dir, "$(species_name)_$(prop)_test.png")
        savefig(p, output_file)
        println("  Saved to $output_file")
    end
    
    # Generate combined plot
    println("  Generating combined plot...")
    p = plot_all_properties(result)
    output_file = joinpath(output_dir, "$(species_name)_all_properties_test.png")
    savefig(p, output_file)
    println("  Saved to $output_file")
end

println("Done!")