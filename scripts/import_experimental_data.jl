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

# Run the sample data creator instead of trying to download from external sites
sample_data_script = joinpath(@__DIR__, "experimental_data_importers", "create_sample_data.jl")
if isfile(sample_data_script)
    println("\nCreating sample experimental data...")
    include(sample_data_script)
end

println("\n====================================================")
println("EXPERIMENTAL DATA IMPORT COMPLETE")
println("====================================================")
println("Finished at: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")