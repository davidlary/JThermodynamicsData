#!/usr/bin/env julia

"""
JThermodynamicsData Command Line Interface

This script provides a command-line interface for JThermodynamicsData.
"""

using Pkg
Pkg.activate(@__DIR__)

using ArgParse
using JThermodynamicsData
using YAML
using DataFrames
using DuckDB
using CSV
using Plots
using Printf

function parse_command_line()
    s = ArgParseSettings(
        description="JThermodynamicsData - A Julia package for thermodynamic calculations",
        version="0.1.0",
        add_version=true
    )
    
    @add_arg_table s begin
        "--config", "-c"
            help = "Path to configuration file"
            default = joinpath(@__DIR__, "config", "settings.yaml")
        
        "command"
            help = """Command to run:
                    - query: Query thermodynamic properties for a species
                    - list: List available species
                    - update: Update thermodynamic data
                    - calculate: Calculate properties using theoretical methods
                    - plot: Generate plots for thermodynamic properties
                    - convert: Convert between data formats"""
            required = true
        
        "args"
            help = "Arguments for the command"
            nargs = '*'
    end
    
    return parse_args(s)
end

function run_query_command(args, config_path)
    if length(args) < 1
        println("Error: Missing species name")
        println("Usage: jthermodynamics.jl query <species> [temperature] [property]")
        return 1
    end
    
    species = args[1]
    
    # Default temperature is 298.15 K
    temperature = 298.15
    if length(args) >= 2
        temperature = parse(Float64, args[2])
    end
    
    # Default property is "all"
    property = "all"
    if length(args) >= 3
        property = args[3]
    end
    
    # Load configuration
    config = load_config(config_path)
    
    # Initialize database
    conn = init_database(config)
    
    try
        # Query the species
        result = query_properties(conn, species, temperature, config)
        
        # Close database connection
        close_database(conn)
        
        # Display results
        println("\nThermodynamic properties for $(species) at $(temperature) K:")
        
        if isempty(result)
            println("No data found for species: $(species)")
            return 1
        end
        
        if property == "all"
            # Display all properties
            for (prop_name, prop_data) in result["properties"]
                value = prop_data["value"]
                uncertainty = prop_data["uncertainty"]
                units = prop_data["units"]
                
                println("  $(prop_name): $(value) ± $(uncertainty) $(units)")
            end
            
            # Display source information
            println("\nData sources used:")
            for source in result["sources"]
                println("  - $(source)")
            end
        else
            # Display specific property
            if haskey(result["properties"], property)
                prop_data = result["properties"][property]
                value = prop_data["value"]
                uncertainty = prop_data["uncertainty"]
                units = prop_data["units"]
                
                println("  $(property): $(value) ± $(uncertainty) $(units)")
            else
                println("Property not found: $(property)")
                println("Available properties: $(keys(result["properties"]))")
                return 1
            end
        end
        
        return 0
    catch e
        # Close database connection
        close_database(conn)
        
        println("Error: $(e)")
        return 1
    end
end

function run_list_command(args, config_path)
    # Load configuration
    config = load_config(config_path)
    
    # Initialize database
    conn = init_database(config)
    
    try
        # Filter option
        filter = ""
        if length(args) >= 1
            filter = args[1]
        end
        
        # Source option
        source = ""
        if length(args) >= 2
            source = args[2]
        end
        
        # Query available species
        if isempty(source)
            # List all species
            query = """
            SELECT s.name, s.formula, s.cas_number, COUNT(td.id) AS data_count
            FROM species s
            LEFT JOIN thermodynamic_data td ON s.id = td.species_id
            """
            
            if !isempty(filter)
                query *= " WHERE s.name LIKE '%$(filter)%' OR s.formula LIKE '%$(filter)%'"
            end
            
            query *= " GROUP BY s.name, s.formula, s.cas_number ORDER BY s.name"
            
            result = DuckDB.execute(conn, query)
            species_df = DataFrame(result)
            
            # Display results
            println("\nAvailable species:")
            if size(species_df, 1) > 0
                for row in eachrow(species_df)
                    println("  - $(row.name) ($(row.formula)), CAS: $(row.cas_number), Data entries: $(row.data_count)")
                end
                println("\nTotal: $(size(species_df, 1)) species")
            else
                println("  No species found")
            end
        else
            # List species for a specific source
            query = """
            SELECT s.name, s.formula, s.cas_number
            FROM species s
            JOIN thermodynamic_data td ON s.id = td.species_id
            WHERE td.data_source = ?
            """
            
            if !isempty(filter)
                query *= " AND (s.name LIKE '%$(filter)%' OR s.formula LIKE '%$(filter)%')"
            end
            
            query *= " ORDER BY s.name"
            
            result = DuckDB.execute(conn, query, [source])
            species_df = DataFrame(result)
            
            # Display results
            println("\nSpecies from source $(source):")
            if size(species_df, 1) > 0
                for row in eachrow(species_df)
                    println("  - $(row.name) ($(row.formula)), CAS: $(row.cas_number)")
                end
                println("\nTotal: $(size(species_df, 1)) species")
            else
                println("  No species found for source: $(source)")
            end
        end
        
        # List available sources
        query = """
        SELECT name, priority, last_updated
        FROM data_sources
        ORDER BY priority DESC
        """
        
        result = DuckDB.execute(conn, query)
        sources_df = DataFrame(result)
        
        println("\nAvailable data sources:")
        if size(sources_df, 1) > 0
            for row in eachrow(sources_df)
                last_updated = ismissing(row.last_updated) ? "Never" : row.last_updated
                println("  - $(row.name) (Priority: $(row.priority)), Last updated: $(last_updated)")
            end
        else
            println("  No data sources found")
        end
        
        # Close database connection
        close_database(conn)
        
        return 0
    catch e
        # Close database connection
        close_database(conn)
        
        println("Error: $(e)")
        return 1
    end
end

function run_update_command(args, config_path)
    # Load configuration
    config = load_config(config_path)
    
    # Specific source to update
    source_name = ""
    if length(args) >= 1
        source_name = args[1]
    end
    
    # Force update flag
    force = false
    if length(args) >= 2 && args[2] == "force"
        force = true
    end
    
    # Call the update_database.jl script
    update_script = joinpath(@__DIR__, "scripts", "update_database.jl")
    
    if isempty(source_name)
        # Update all sources
        cmd = `julia $(update_script)`
        run(cmd)
    else
        # Update specific source
        cmd = `julia $(update_script) $(source_name) $(force ? "force" : "")`
        run(cmd)
    end
    
    return 0
end

function run_calculate_command(args, config_path)
    if length(args) < 1
        println("Error: Missing species formula")
        println("Usage: jthermodynamics.jl calculate <formula> [temperature] [method]")
        return 1
    end
    
    formula = args[1]
    
    # Default temperature is 298.15 K
    temperature = 298.15
    if length(args) >= 2
        temperature = parse(Float64, args[2])
    end
    
    # Default method is "all"
    method = "all"
    if length(args) >= 3
        method = args[3]
    end
    
    # Load configuration
    config = load_config(config_path)
    
    try
        # Calculate properties using theoretical methods
        println("\nCalculating properties for $(formula) at $(temperature) K using $(method)...")
        
        if method == "all"
            # Try all methods
            methods = ["group_contribution", "statistical_thermodynamics", "machine_learning"]
            
            for m in methods
                try
                    if m == "group_contribution"
                        result = estimate_properties_group_contribution(formula, temperature)
                        display_calculated_properties(result, m)
                    elseif m == "statistical_thermodynamics"
                        # Estimate molecular properties first
                        molecular_data = estimate_molecular_properties(formula)
                        result = estimate_properties_statistical_thermodynamics(formula, molecular_data, temperature)
                        display_calculated_properties(result, m)
                    elseif m == "machine_learning"
                        result = estimate_properties_machine_learning(formula, temperature, config)
                        display_calculated_properties(result, m)
                    end
                catch e
                    println("  Method $(m) failed: $(e)")
                end
            end
        else
            # Try specific method
            if method == "group_contribution"
                result = estimate_properties_group_contribution(formula, temperature)
                display_calculated_properties(result, method)
            elseif method == "statistical_thermodynamics"
                # Estimate molecular properties first
                molecular_data = estimate_molecular_properties(formula)
                result = estimate_properties_statistical_thermodynamics(formula, molecular_data, temperature)
                display_calculated_properties(result, method)
            elseif method == "machine_learning"
                result = estimate_properties_machine_learning(formula, temperature, config)
                display_calculated_properties(result, method)
            elseif method == "joback"
                # For Joback method, we need SMILES
                println("  Joback method requires SMILES notation. Using formula as approximation.")
                result = joback_method(formula, temperature)
                display_calculated_properties(result, method)
            else
                println("Unknown method: $(method)")
                println("Available methods: group_contribution, statistical_thermodynamics, machine_learning, joback")
                return 1
            end
        end
        
        return 0
    catch e
        println("Error: $(e)")
        return 1
    end
end

function display_calculated_properties(result, method)
    println("\nResults from $(method):")
    println("  Cp: $(result["Cp"]) ± $(result["Cp_uncertainty"]) J/mol/K")
    println("  H: $(result["H"]) ± $(result["H_uncertainty"]) kJ/mol")
    println("  S: $(result["S"]) ± $(result["S_uncertainty"]) J/mol/K")
    println("  G: $(result["G"]) ± $(result["G_uncertainty"]) kJ/mol")
end

function run_plot_command(args, config_path)
    if length(args) < 1
        println("Error: Missing species name")
        println("Usage: jthermodynamics.jl plot <species> [property] [temp_min] [temp_max]")
        return 1
    end
    
    species = args[1]
    
    # Default property is "all"
    property = "all"
    if length(args) >= 2
        property = args[2]
    end
    
    # Default temperature range
    temp_min = 100.0
    temp_max = 3000.0
    
    if length(args) >= 4
        temp_min = parse(Float64, args[3])
        temp_max = parse(Float64, args[4])
    end
    
    # Load configuration
    config = load_config(config_path)
    
    # Initialize database
    conn = init_database(config)
    
    try
        # Generate temperature points (logarithmic scale)
        n_points = 100
        temperatures = exp.(range(log(temp_min), log(temp_max), length=n_points))
        
        # Calculate properties
        println("\nCalculating properties for $(species) over temperature range $(temp_min)-$(temp_max) K...")
        
        property_data = Dict(
            "species_name" => species,
            "temperatures" => temperatures,
            "properties" => Dict()
        )
        
        # Calculate properties at each temperature
        for temp in temperatures
            result = query_properties(conn, species, temp, config)
            
            if isempty(result)
                println("No data found for species: $(species)")
                return 1
            end
            
            # Store values for each property
            for (prop_name, prop_data) in result["properties"]
                if !haskey(property_data["properties"], prop_name)
                    property_data["properties"][prop_name] = Dict(
                        "values" => Float64[],
                        "uncertainties" => Float64[],
                        "units" => prop_data["units"]
                    )
                end
                
                push!(property_data["properties"][prop_name]["values"], prop_data["value"])
                push!(property_data["properties"][prop_name]["uncertainties"], prop_data["uncertainty"])
            end
        end
        
        # Close database connection
        close_database(conn)
        
        # Generate plots
        if property == "all"
            # Create a 2x2 plot with all properties
            p_cp = plot_property(property_data, "Cp", temperatures, true)
            p_h = plot_property(property_data, "H", temperatures, true)
            p_s = plot_property(property_data, "S", temperatures, true)
            p_g = plot_property(property_data, "G", temperatures, true)
            
            p = plot(p_cp, p_h, p_s, p_g, layout=(2, 2), size=(1000, 800),
                    title="Thermodynamic Properties for $(species)")
            
            # Save plot
            output_file = "$(species)_properties.png"
            savefig(p, output_file)
            println("Plot saved to: $(output_file)")
        else
            # Create a single plot for the specified property
            p = plot_property(property_data, property, temperatures, false)
            
            # Save plot
            output_file = "$(species)_$(property).png"
            savefig(p, output_file)
            println("Plot saved to: $(output_file)")
        end
        
        return 0
    catch e
        # Close database connection
        close_database(conn)
        
        println("Error: $(e)")
        return 1
    end
end

function plot_property(data, property, temperatures, as_subplot)
    if !haskey(data["properties"], property)
        error("Property not found: $(property)")
    end
    
    prop_data = data["properties"][property]
    values = prop_data["values"]
    uncertainties = prop_data["uncertainties"]
    units = prop_data["units"]
    
    title = as_subplot ? property : "$(property) for $(data["species_name"])"
    
    # Create plot with uncertainty ribbon
    p = plot(temperatures, values,
            xscale=:log10,
            xlabel="Temperature (K)",
            ylabel="$(property) ($(units))",
            title=title,
            ribbon=uncertainties,
            fillalpha=0.3,
            legend=false,
            lw=2,
            grid=true,
            framestyle=:box)
    
    return p
end

function run_convert_command(args, config_path)
    if length(args) < 3
        println("Error: Missing arguments")
        println("Usage: jthermodynamics.jl convert <input_file> <output_file> <format>")
        return 1
    end
    
    input_file = args[1]
    output_file = args[2]
    format = args[3]
    
    # Load configuration
    config = load_config(config_path)
    
    try
        # Convert data format
        println("\nConverting $(input_file) to format: $(format)...")
        
        if !isfile(input_file)
            println("Input file not found: $(input_file)")
            return 1
        end
        
        if format == "nasa7"
            # Convert to NASA 7-coefficient format
            convert_to_nasa7(input_file, output_file)
        elseif format == "nasa9"
            # Convert to NASA 9-coefficient format
            convert_to_nasa9(input_file, output_file)
        elseif format == "shomate"
            # Convert to Shomate format
            convert_to_shomate(input_file, output_file)
        elseif format == "csv"
            # Convert to CSV format
            convert_to_csv(input_file, output_file)
        else
            println("Unsupported output format: $(format)")
            println("Supported formats: nasa7, nasa9, shomate, csv")
            return 1
        end
        
        println("Conversion complete. Output saved to: $(output_file)")
        
        return 0
    catch e
        println("Error: $(e)")
        return 1
    end
end

function convert_to_nasa7(input_file, output_file)
    # Implement conversion to NASA 7-coefficient format
    # This is a placeholder for the actual implementation
    error("Conversion to NASA 7-coefficient format not implemented yet")
end

function convert_to_nasa9(input_file, output_file)
    # Implement conversion to NASA 9-coefficient format
    # This is a placeholder for the actual implementation
    error("Conversion to NASA 9-coefficient format not implemented yet")
end

function convert_to_shomate(input_file, output_file)
    # Implement conversion to Shomate format
    # This is a placeholder for the actual implementation
    error("Conversion to Shomate format not implemented yet")
end

function convert_to_csv(input_file, output_file)
    # Implement conversion to CSV format
    # This is a placeholder for the actual implementation
    error("Conversion to CSV format not implemented yet")
end

function main()
    args = parse_command_line()
    
    command = args["command"]
    command_args = args["args"]
    config_path = args["config"]
    
    # Run the appropriate command
    exit_code = 0
    
    if command == "query"
        exit_code = run_query_command(command_args, config_path)
    elseif command == "list"
        exit_code = run_list_command(command_args, config_path)
    elseif command == "update"
        exit_code = run_update_command(command_args, config_path)
    elseif command == "calculate"
        exit_code = run_calculate_command(command_args, config_path)
    elseif command == "plot"
        exit_code = run_plot_command(command_args, config_path)
    elseif command == "convert"
        exit_code = run_convert_command(command_args, config_path)
    else
        println("Unknown command: $(command)")
        println("Available commands: query, list, update, calculate, plot, convert")
        exit_code = 1
    end
    
    return exit_code
end

# Run the main function and exit with the appropriate code
exit(main())