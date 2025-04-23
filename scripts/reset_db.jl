#!/usr/bin/env julia

# Reset Database Script
# This script creates a fresh database by removing the existing one

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using DuckDB
using YAML

# Load settings
config_path = joinpath(dirname(dirname(@__FILE__)), "config", "settings.yaml")
config = YAML.load_file(config_path; dicttype=Dict{String,Any})

# Get database path
db_path = config["general"]["database_path"]

# Remove existing database if it exists
if isfile(db_path)
    println("Removing existing database: $(db_path)")
    rm(db_path)
    println("Database removed successfully.")
else
    println("No existing database found at: $(db_path)")
end

# Ensure directory exists
mkpath(dirname(db_path))
println("Ready to create a new database.")