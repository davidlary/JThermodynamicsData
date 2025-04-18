#!/usr/bin/env julia

"""
JThermodynamicsData Setup Script

This script installs all required dependencies for JThermodynamicsData.
Run this before using the main program.
"""

using Pkg
println("Activating project at $(dirname(@__FILE__))")
Pkg.activate(dirname(@__FILE__))

# Install required packages
println("Installing required packages...")
packages = [
    "YAML",
    "DuckDB",
    "HTTP",
    "CSV",
    "JSON",
    "DataFrames",
    "Measurements",
    "MonteCarloMeasurements",
    "Interpolations", 
    "Statistics",
    "Plots",
    "Unitful",
    "UnitfulMoles",
    "LightXML",
    "StaticArrays",
    "SpecialFunctions",
    "ArgParse"
]

for pkg in packages
    println("Installing $pkg...")
    Pkg.add(pkg)
end

# Create required directories
println("Creating required directories...")
mkpath(joinpath(dirname(@__FILE__), "data", "cache"))
mkpath(joinpath(dirname(@__FILE__), "output", "tables"))
mkpath(joinpath(dirname(@__FILE__), "plots"))

println("\nSetup complete! You can now run the main program or CLI:")
println("julia run_thermodynamics.jl")
println("./jthermodynamicsdata.jl --help")