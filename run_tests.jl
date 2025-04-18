#!/usr/bin/env julia

"""
Run all tests for JThermodynamicsData package

This script runs the test_parsers.jl script to download sample data and test all parsers.
"""

using Pkg
Pkg.activate(dirname(@__FILE__))

# Run the parser tests
include(joinpath(dirname(@__FILE__), "scripts", "test_parsers.jl"))

println("\nAll tests completed!")