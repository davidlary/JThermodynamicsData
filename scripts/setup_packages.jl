#!/usr/bin/env julia

"""
Setup Packages Script

This script ensures all necessary packages are installed for JThermodynamicsData.
Run this script before running any other scripts in the package.

Usage:
    julia scripts/setup_packages.jl
"""

using Pkg

# Main project packages
Pkg.activate(".")
println("Installing main project packages...")

packages = [
    "CSV",
    "Dates",
    "DocStringExtensions",
    "DuckDB", 
    "JLD2",
    "HTTP",
    "LsqFit",
    "Printf",
    "YAML",
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
    "LinearAlgebra",
    "DelimitedFiles",
    "Optim"
]

for pkg in packages
    println("Adding $pkg...")
    try
        Pkg.add(pkg)
    catch e
        if isa(e, Pkg.Types.PkgError)
            println("  Package $pkg is in standard library or already installed.")
        else
            println("  Error installing $pkg: $(e)")
        end
    end
end

# Scripts directory packages
Pkg.activate("./scripts")
println("\nInstalling script-specific packages...")

script_packages = [
    "DelimitedFiles",
    "Statistics",
    "Plots",
    "Printf",
    "Dates",
    "YAML",
    "LinearAlgebra",
    "HTTP",
    "SpecialFunctions"
]

for pkg in script_packages
    println("Adding $pkg to scripts environment...")
    try
        Pkg.add(pkg)
    catch e
        if isa(e, Pkg.Types.PkgError)
            println("  Package $pkg is in standard library or already installed.")
        else
            println("  Error installing $pkg: $(e)")
        end
    end
end

println("\nAll packages installed successfully!")
println("You can now run the scripts without package errors.")