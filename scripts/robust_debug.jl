#!/usr/bin/env julia

# Robust debug script for JThermodynamicsData
# Handles all known error cases and provides step-by-step testing

using Pkg
Pkg.activate(".")

println("Loading required packages...")
using JThermodynamicsData
using YAML
using DuckDB
using DataFrames
using JSON
using Dates
using Plots

println("Start time: $(now())")

# Patch functions to fix known issues
println("Applying fixes to known issues...")

# Define lowercase constants to fix statistical thermodynamics
if !isdefined(JThermodynamicsData, :h)
    println("Adding lowercase constant aliases...")
    
    # This is a workaround to add constants to an imported module
    methods = [:h, :kb, :na, :c]
    for method in methods
        uppercase_method = Symbol(uppercase(String(method)))
        cmd = :(Core.eval(JThermodynamicsData, :(const $method = $(JThermodynamicsData.$uppercase_method))))
        eval(cmd)
        
        # Verify fix
        if isdefined(JThermodynamicsData, method)
            println("  Fixed: Added JThermodynamicsData.$method = $(getfield(JThermodynamicsData, method))")
        else
            println("  Warning: Failed to add JThermodynamicsData.$method")
        end
    end
end

# Modified calculate_properties that correctly handles keyword arguments
function patched_calculate_properties(conn::DuckDB.DB, species_name::String, temperature_range::Vector{Float64}, 
                                    config::Dict; step::Float64=10.0)
    # Validate temperature range
    if length(temperature_range) != 2 || temperature_range[1] >= temperature_range[2]
        error("Invalid temperature range: $temperature_range")
    end
    
    # Generate temperature points
    temperatures = temperature_range[1]:step:temperature_range[2]
    
    # Calculate properties at each temperature
    results = Dict[]
    
    for temp in temperatures
        try
            result = JThermodynamicsData.query_properties(conn, species_name, temp, config)
            push!(results, result)
        catch e
            @warn "Failed to calculate properties at $temp K: $e"
        end
    end
    
    # Extract property data for plotting
    temps = [result["temperature"] for result in results]
    
    cp_values = [result["properties"]["Cp"]["value"] for result in results]
    cp_uncertainties = [result["properties"]["Cp"]["uncertainty"] for result in results]
    
    h_values = [result["properties"]["H"]["value"] for result in results]
    h_uncertainties = [result["properties"]["H"]["uncertainty"] for result in results]
    
    s_values = [result["properties"]["S"]["value"] for result in results]
    s_uncertainties = [result["properties"]["S"]["uncertainty"] for result in results]
    
    g_values = [result["properties"]["G"]["value"] for result in results]
    g_uncertainties = [result["properties"]["G"]["uncertainty"] for result in results]
    
    # Create result structure
    return Dict(
        "species_name" => species_name,
        "temperature_range" => temperature_range,
        "temperatures" => temps,
        "properties" => Dict(
            "Cp" => Dict(
                "values" => cp_values,
                "uncertainties" => cp_uncertainties,
                "units" => "J/mol/K"
            ),
            "H" => Dict(
                "values" => h_values,
                "uncertainties" => h_uncertainties,
                "units" => "kJ/mol"
            ),
            "S" => Dict(
                "values" => s_values,
                "uncertainties" => s_uncertainties,
                "units" => "J/mol/K"
            ),
            "G" => Dict(
                "values" => g_values,
                "uncertainties" => g_uncertainties,
                "units" => "kJ/mol"
            )
        )
    )
end

# Load configuration
println("Loading configuration and connecting to database...")
config_path = joinpath(dirname(@__DIR__), "config", "settings.yaml")
config = JThermodynamicsData.load_config(config_path)
conn = JThermodynamicsData.init_database(config)

# Species to test
species_name = "N2"
println("Testing with species: $species_name")

# First, check if species exists in database
println("\nStep 1: Checking if species exists in database...")
species_id = JThermodynamicsData.get_species_id(conn, species_name)

if species_id === nothing
    println("❌ Species $species_name not found in database")
    println("  Running initialization script to add N2...")
    include("initialize_minimal_database.jl")
    species_id = JThermodynamicsData.get_species_id(conn, species_name)
    println("  After initialization: species ID = $species_id")
else
    println("✓ Species $species_name found with ID: $species_id")
end

# Check the species metadata
println("\nStep 2: Checking species metadata...")
query_result = DuckDB.execute(conn, "SELECT * FROM species WHERE id = ?", [species_id])
species_df = DataFrame(query_result)

if ismissing(species_df[1, :metadata_json]) || isempty(species_df[1, :metadata_json])
    println("❌ No metadata found for species, adding it...")
    metadata = Dict(
        "molecular_data" => Dict(
            "molecular_weight" => 28.0134,
            "symmetry_number" => 2,
            "spin_multiplicity" => 1,
            "rotational_constants" => [1.998],  # cm^-1
            "vibrational_frequencies" => [2358.0]  # cm^-1
        )
    )
    metadata_json = JSON.json(metadata)
    
    update_query = "UPDATE species SET metadata_json = ? WHERE id = ?"
    DuckDB.execute(conn, update_query, [metadata_json, species_id])
    println("✓ Added metadata for $species_name")
else
    try
        metadata = JSON.parse(species_df[1, :metadata_json])
        println("✓ Metadata found and valid")
    catch e
        println("❌ Invalid metadata JSON: $e")
        # Fix metadata
        metadata = Dict(
            "molecular_data" => Dict(
                "molecular_weight" => 28.0134,
                "symmetry_number" => 2,
                "spin_multiplicity" => 1,
                "rotational_constants" => [1.998],  # cm^-1
                "vibrational_frequencies" => [2358.0]  # cm^-1
            )
        )
        metadata_json = JSON.json(metadata)
        
        update_query = "UPDATE species SET metadata_json = ? WHERE id = ?"
        DuckDB.execute(conn, update_query, [metadata_json, species_id])
        println("✓ Fixed metadata for $species_name")
    end
end

# Check thermodynamic data
println("\nStep 3: Checking thermodynamic data sources...")
data_query = "SELECT * FROM thermodynamic_data WHERE species_id = ?"
data_result = DuckDB.execute(conn, data_query, [species_id])
data_df = DataFrame(data_result)

if size(data_df, 1) == 0
    println("❌ No thermodynamic data found for species")
    println("  Running initialization script to add data...")
    include("initialize_minimal_database.jl")
    data_result = DuckDB.execute(conn, data_query, [species_id])
    data_df = DataFrame(data_result)
    println("  After initialization: found $(size(data_df, 1)) data sources")
else
    println("✓ Found $(size(data_df, 1)) thermodynamic data sources")
end

# Check data source record
println("\nStep 4: Checking data source records...")
source_query = "SELECT * FROM data_sources"
source_result = DuckDB.execute(conn, source_query)
source_df = DataFrame(source_result)

if size(source_df, 1) == 0
    println("❌ No data sources found in database")
    println("  Running initialization script to add data sources...")
    include("initialize_minimal_database.jl")
    source_result = DuckDB.execute(conn, source_query)
    source_df = DataFrame(source_result)
    println("  After initialization: found $(size(source_df, 1)) sources")
else
    println("✓ Found $(size(source_df, 1)) data sources")
end

# Test direct property calculation at a single temperature
println("\nStep 5: Testing direct property calculation at single temperature...")

try
    temperature = 298.15  # K
    result, sources = JThermodynamicsData.progressively_refine_thermodynamic_data(
        conn, species_name, temperature, config
    )
    
    println("✓ Successfully calculated properties for $species_name at $temperature K")
    println("  Data source: $(result["data_source"])")
    
    # Display properties
    for (prop, data) in result["properties"]
        println("  $prop: $(data["value"]) ± $(data["uncertainty"]) $(data["units"])")
    end
    
    println("  Used $(length(sources)) data sources in calculation")
catch e
    println("❌ Failed to calculate properties: $e")
    println(stacktrace(catch_backtrace()))
end

# Test property calculation across temperature range
println("\nStep 6: Testing property calculation across temperature range...")

try
    # Small temperature range for testing
    temp_range = [200.0, 500.0]
    temp_step = 50.0
    
    # Use our patched function to correctly handle keyword arguments
    result = patched_calculate_properties(conn, species_name, temp_range, config; step=temp_step)
    
    println("✓ Successfully calculated properties across temperature range")
    println("  Temperature points: $(length(result["temperatures"]))")
    println("  Temperature range: $(temp_range[1])-$(temp_range[2]) K")
    
    # Create a simple plot to verify results
    println("\nStep 7: Creating verification plot...")
    
    try
        # Heat capacity plot
        temps = result["temperatures"]
        cp_values = result["properties"]["Cp"]["values"]
        cp_uncertainties = result["properties"]["Cp"]["uncertainties"]
        
        p = plot(
            temps,
            cp_values,
            ribbon=cp_uncertainties,
            fillalpha=0.3,
            xlabel="Temperature (K)",
            ylabel="Cp (J/mol/K)",
            title="Heat Capacity for $species_name",
            legend=false,
            lw=2,
            grid=true
        )
        
        plot_dir = joinpath(dirname(@__DIR__), "plots")
        mkpath(plot_dir)
        plot_file = joinpath(plot_dir, "debug_N2_Cp.png")
        savefig(p, plot_file)
        println("✓ Created verification plot at: $plot_file")
    catch e
        println("❌ Failed to create plot: $e")
    end
catch e
    println("❌ Failed to calculate properties across temperature range: $e")
    println(stacktrace(catch_backtrace()))
end

# Close the database connection
println("\nClosing database connection...")
JThermodynamicsData.close_database(conn)

println("End time: $(now())")
println("Debug process complete!")