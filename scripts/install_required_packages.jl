#!/usr/bin/env julia

"""
Install Required Packages

This script installs all required packages for the JThermodynamicsData workflow.
"""

using Pkg

println("Installing required packages...")

# Activate the project
Pkg.activate(dirname(dirname(@__FILE__)))

# Add required packages for data download and processing
packages = [
    "HTTP",
    "JSON",
    "YAML",
    "ZipFile",
    "Gumbo",
    "Cascadia",
    "Tar",
    "CodecZlib",
    "LazyArtifacts",
    "DataFrames",
    "DuckDB",
    "CSV",
    "Plots",
    "Colors",
    "Measurements",
    "MonteCarloMeasurements",
    "Interpolations",
    "Statistics",
    "Unitful",
    "UnitfulMoles",
    "LightXML",
    "StaticArrays",
    "SpecialFunctions"
]

for pkg in packages
    println("Adding package: $pkg")
    Pkg.add(pkg)
end

println("All required packages installed.")
println("\nNow you can run the full workflow:")
println("julia run_full_workflow.jl")