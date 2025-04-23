#!/usr/bin/env julia

# Hierarchical Thermodynamics Demonstration Runner
# This script runs the hierarchical thermodynamic calculator demonstration

println("JThermodynamicsData - Hierarchical Calculator Demonstration")
println("===========================================================")

# Activate the project
using Pkg
Pkg.activate(@__DIR__)

try
    # Make sure the required packages are installed
    Pkg.instantiate()
    
    # Run the demonstration script
    include(joinpath(@__DIR__, "examples/hierarchical_demo.jl"))
    
    # Run the main function in the demonstration
    main()
    
    println("\nDemonstration completed successfully!")
    println("Check the 'plots' directory for visualization results.")
catch e
    println("\nError running the demonstration: ", e)
    
    # Detailed error information
    println("\nError details:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end