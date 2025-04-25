#!/usr/bin/env julia

"""
Import Experimental Data

This script imports data from various experimental sources:
- ATCT (Active Thermochemical Tables)
- JANAF Thermochemical Tables
- NASA CEA
- NIST Webbook
- ThermoML
- TDE (ThermoData Engine)
"""

using Pkg
Pkg.activate(@__DIR__)
using Dates

println("\n====================================================")
println("IMPORTING EXPERIMENTAL THERMODYNAMIC DATA")
println("====================================================")
println("Starting at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")

# Run ATCT importer
atct_script = joinpath(@__DIR__, "experimental_data_importers", "atct_importer.jl")
if isfile(atct_script)
    println("\nRunning ATCT importer...")
    include(atct_script)
end

# Run JANAF importer
janaf_script = joinpath(@__DIR__, "experimental_data_importers", "janaf_importer.jl")
if isfile(janaf_script)
    println("\nRunning JANAF importer...")
    include(janaf_script)
end

println("\n====================================================")
println("EXPERIMENTAL DATA IMPORT COMPLETE")
println("====================================================")
println("Finished at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")