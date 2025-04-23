#!/usr/bin/env julia

# Direct Calculator Runner
# This script runs the direct calculator demonstration with hierarchical data fetching,
# logarithmic temperature scales, and uncertainty visualization.

println("JThermodynamicsData - Direct Calculator Demonstration")
println("=====================================================")
println("Running with: Log temperature scales & uncertainty visualization")

# Activate the project
using Pkg
Pkg.activate(@__DIR__)

try
    # Make sure the required packages are installed
    Pkg.instantiate()
    
    # Run the demonstration script
    include(joinpath(@__DIR__, "examples/direct_calculator.jl"))
    
    println("\nDemonstration completed successfully!")
    println("\nOutputs:")
    println("  - Plots: ./plots/[species]_all_properties.png")
    println("  - Data tables: ./output/tables/[species]_hierarchical.csv")
    println("  - NASA-7 polynomials: ./output/tables/[species]_nasa7.txt")
    println("  - Source documentation: ./output/docs/[species]_sources.md")
    println("  - Summary: ./output/tables/all_species_summary.csv")
catch e
    println("\nError running the demonstration: ", e)
    
    # Detailed error information
    println("\nError details:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end