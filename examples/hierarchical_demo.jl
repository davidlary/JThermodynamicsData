#!/usr/bin/env julia

# Hierarchical Thermodynamic Calculator Demonstration
# This script demonstrates how the system selects thermodynamic data from
# different sources based on configured priority.

using Pkg
# Activate the current project
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using YAML
using DataFrames
using Plots
using Colors
using Printf

# Include sample data for demonstration purposes
include(joinpath(dirname(dirname(@__FILE__)), "data/external/demo/sample_data.jl"))

function initialize_demo_data(db)
    # Clear any existing data for demo purposes
    println("Initializing demo data...")
    conn = JThermodynamicsData.get_connection(db)
    DBInterface.execute(conn, "DELETE FROM thermodynamic_data")
    DBInterface.execute(conn, "DELETE FROM species")
    DBInterface.execute(conn, "DELETE FROM data_sources")
    
    # Add our sample data sources
    sources = [
        ("theoretical", "Theoretical calculations", 1, 3.0),
        ("gri-mech", "GRI-MECH 3.0", 3, 3.5),
        ("burcat", "Burcat database", 7, 4.2),
        ("janaf", "JANAF tables", 8, 4.5),
        ("atct", "Active Thermochemical Tables", 9, 5.0)
    ]
    
    for (i, (name, description, priority, reliability)) in enumerate(sources)
        # Create metadata with reliability score
        metadata = Dict(
            "format" => "nasa7",
            "enabled" => true,
            "reliability_score" => reliability,
            "refresh_interval_days" => 90
        )
        
        metadata_json = JSON.json(metadata)
        
        DBInterface.execute(conn, """
            INSERT INTO data_sources (id, name, description, priority, metadata_json)
            VALUES (?, ?, ?, ?, ?)
        """, [i, name, description, priority, metadata_json])
    end
    
    # Add the species from the sample data
    for source_name in keys(SAMPLE_DATA)
        # Add each species for this source
        for species_name in keys(SAMPLE_DATA[source_name])
            species_data = SAMPLE_DATA[source_name][species_name]
            
            # Check if species exists, if not add it
            result = DBInterface.execute(conn, "SELECT id FROM species WHERE name = ?", [species_name])
            if size(DataFrame(result), 1) == 0
                DBInterface.execute(conn, """
                    INSERT INTO species (id, name, formula, cas_number, molecular_weight)
                    VALUES (nextval('species_id_seq'), ?, ?, ?, 0.0)
                """, [species_name, species_data["formula"], species_data["cas"]])
            end
            
            # Get the species ID
            result = DBInterface.execute(conn, "SELECT id FROM species WHERE name = ?", [species_name])
            species_id = DataFrame(result)[1, :id]
            
            # We need to add data for both temperature ranges
            temp_ranges = species_data["temperature_ranges"]
            coeffs = species_data["coefficients"]
            
            # Convert coefficients to JSON data
            data_json = JSON.json(Dict(
                "low_temp" => Dict(
                    "range" => temp_ranges[1],
                    "coefficients" => coeffs[1]
                ),
                "high_temp" => Dict(
                    "range" => temp_ranges[2],
                    "coefficients" => coeffs[2]
                )
            ))
            
            # Convert uncertainty to JSON
            uncertainty_json = JSON.json(Dict(
                "type" => "constant",
                "value" => get(species_data, "reliability", 0.1)
            ))
            
            # Add thermodynamic data
            DBInterface.execute(conn, """
                INSERT INTO thermodynamic_data 
                (id, species_id, data_source, polynomial_type, temperature_min, temperature_max, 
                 reliability_score, data_json, uncertainty_json)
                VALUES (nextval('thermodynamic_data_id_seq'), ?, ?, 'nasa7', ?, ?, ?, ?, ?)
            """, [
                species_id, source_name, 
                temp_ranges[1][1], temp_ranges[2][2],
                species_data["reliability"],
                data_json, uncertainty_json
            ])
        end
    end
    
    println("Demo data loaded successfully!")
end

function test_hierarchical_selection(db, species, temp=298.15)
    println("\n===== Testing hierarchical selection for $(species) at $(temp)K =====")
    
    # Get available sources for this species
    conn = JThermodynamicsData.get_connection(db)
    query = """
        SELECT ds.name as source, ds.priority, 
               CAST(json_extract(ds.metadata_json, '$.reliability_score') AS DOUBLE) as reliability_score
        FROM thermodynamic_data td
        JOIN data_sources ds ON td.data_source = ds.name
        JOIN species sp ON td.species_id = sp.id
        WHERE sp.name = ?
        ORDER BY ds.priority
    """
    result = DBInterface.execute(conn, query, [species])
    sources_df = DataFrame(result)
    
    if size(sources_df, 1) == 0
        println("No data found for species: $(species)")
        return
    end
    
    println("Available sources:")
    for row in eachrow(sources_df)
        println("  - $(rpad(row.source, 12)) (priority: $(row.priority), reliability: $(row.reliability_score))")
    end
    
    # Calculate properties from each source directly
    source_results = Dict()
    for source_name in sources_df.source
        try
            # Override the hierarchical selection by specifying the source
            poly = JThermodynamicsData.load_polynomial_data(db, species, [temp-1, temp+1], Dict("source" => source_name))
            cp = JThermodynamicsData.calculate_cp(poly, temp)
            h = JThermodynamicsData.calculate_h(poly, temp)
            s = JThermodynamicsData.calculate_s(poly, temp)
            g = JThermodynamicsData.calculate_g(poly, temp)
            
            source_results[source_name] = (cp=cp, h=h, s=s, g=g)
            
            println("\nResults from $(source_name):")
            println("  Cp = $(@sprintf("%.2f", cp)) J/mol/K")
            println("  H  = $(@sprintf("%.2f", h)) kJ/mol")
            println("  S  = $(@sprintf("%.2f", s)) J/mol/K")
            println("  G  = $(@sprintf("%.2f", g)) kJ/mol")
        catch e
            println("Error calculating properties from $(source_name): $(e)")
        end
    end
    
    # Calculate using hierarchical selection (should use highest priority source)
    try
        println("\nResults from hierarchical selection:")
        poly = JThermodynamicsData.load_polynomial_data(db, species, [temp-1, temp+1], Dict("method" => "hierarchical"))
        # Get the selected source
        selected_source = poly.source
        
        cp = JThermodynamicsData.calculate_cp(poly, temp)
        h = JThermodynamicsData.calculate_h(poly, temp)
        s = JThermodynamicsData.calculate_s(poly, temp)
        g = JThermodynamicsData.calculate_g(poly, temp)
        
        println("  Selected source: $(selected_source)")
        println("  Cp = $(@sprintf("%.2f", cp)) J/mol/K")
        println("  H  = $(@sprintf("%.2f", h)) kJ/mol")
        println("  S  = $(@sprintf("%.2f", s)) J/mol/K")
        println("  G  = $(@sprintf("%.2f", g)) kJ/mol")
    catch e
        println("Error calculating properties using hierarchical selection: $(e)")
    end
end

function visualize_results(db, species_list, property="cp", temp_range=[300.0, 1500.0])
    # Create a plot to visualize differences between data sources
    plt = plot(
        title="$(uppercase(property)) Comparison Across Data Sources",
        xlabel="Temperature (K)",
        ylabel=property == "cp" ? "Cp (J/mol/K)" :
               property == "h" ? "H (kJ/mol)" :
               property == "s" ? "S (J/mol/K)" : "G (kJ/mol)",
        legend=:outertopright,
        size=(900, 600),
        dpi=300
    )
    
    # Generate a color palette
    colors = distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    
    # Plot temperature vs property for each species and source
    for (i, species) in enumerate(species_list)
        println("\nVisualizing $(property) for $(species)...")
        
        # Get available sources for this species
        conn = JThermodynamicsData.get_connection(db)
        query = """
            SELECT ds.name as source, ds.priority, 
                   CAST(json_extract(ds.metadata_json, '$.reliability_score') AS DOUBLE) as reliability_score
            FROM thermodynamic_data td
            JOIN data_sources ds ON td.data_source = ds.name
            JOIN species sp ON td.species_id = sp.id
            WHERE sp.name = ?
            ORDER BY ds.priority
        """
        result = DBInterface.execute(conn, query, [species])
        sources_df = DataFrame(result)
        
        if size(sources_df, 1) == 0
            println("No data found for species: $(species)")
            continue
        end
        
        # Calculate hierarchical result
        try
            temps = range(temp_range[1], temp_range[2], length=50)
            poly = JThermodynamicsData.load_polynomial_data(db, species, temp_range, Dict("method" => "hierarchical"))
            
            # Get the property values
            values = if property == "cp"
                [JThermodynamicsData.calculate_cp(poly, t) for t in temps]
            elseif property == "h"
                [JThermodynamicsData.calculate_h(poly, t) for t in temps]
            elseif property == "s"
                [JThermodynamicsData.calculate_s(poly, t) for t in temps]
            else # g
                [JThermodynamicsData.calculate_g(poly, t) for t in temps]
            end
            
            plot!(plt, temps, values, 
                  label="$(species) (hierarchical)", 
                  color=colors[i], 
                  linewidth=3, 
                  linestyle=:solid)
            
            # Add a marker to indicate the selected source
            annotate!(plt, temp_range[2], values[end], 
                      text("$(poly.source)", :right, 8))
        catch e
            println("Error plotting hierarchical data for $(species): $(e)")
        end
        
        # Plot for each source
        for (j, row) in enumerate(eachrow(sources_df))
            try
                source_name = row.source
                temps = range(temp_range[1], temp_range[2], length=50)
                poly = JThermodynamicsData.load_polynomial_data(db, species, temp_range, Dict("source" => source_name))
                
                # Get the property values
                values = if property == "cp"
                    [JThermodynamicsData.calculate_cp(poly, t) for t in temps]
                elseif property == "h"
                    [JThermodynamicsData.calculate_h(poly, t) for t in temps]
                elseif property == "s"
                    [JThermodynamicsData.calculate_s(poly, t) for t in temps]
                else # g
                    [JThermodynamicsData.calculate_g(poly, t) for t in temps]
                end
                
                plot!(plt, temps, values, 
                      label="$(species) ($(source_name))", 
                      color=colors[i], 
                      linewidth=1.5, 
                      linestyle=:dash,
                      alpha=0.7)
            catch e
                println("Error plotting data for $(species) from $(source_name): $(e)")
            end
        end
    end
    
    # Ensure plots directory exists
    plots_dir = joinpath(dirname(dirname(@__FILE__)), "plots")
    if !isdir(plots_dir)
        mkpath(plots_dir)
    end
    
    # Save the plot
    savefig(plt, joinpath(plots_dir, "$(property)_comparison.png"))
    
    println("\nPlot saved to plots/$(property)_comparison.png")
    return plt
end

function main()
    # Initialize the database with demo data
    println("Starting hierarchical calculator demonstration...")
    
    # Load the configuration
    config_path = joinpath(dirname(dirname(@__FILE__)), "config/settings.yaml")
    config = YAML.load_file(config_path; dicttype=Dict{String,Any})
    
    # Initialize database
    db_path = config["general"]["database_path"]
    if !isfile(db_path)
        # Create directory if it doesn't exist
        db_dir = dirname(db_path)
        if !isdir(db_dir)
            mkpath(db_dir)
        end
    end
    
    # Initialize database connection
    db = initialize_database(db_path)
    
    # Load sample data
    initialize_demo_data(db)
    
    # Test individual species
    test_species = ["H2O", "CO2", "N2", "O2"]
    for species in test_species
        test_hierarchical_selection(db, species)
    end
    
    # Visualize the results for different properties
    properties = ["cp", "h", "s", "g"]
    for prop in properties
        visualize_results(db, test_species, prop)
    end
    
    println("\nDemonstration completed successfully!")
end

# Run the main function
main()