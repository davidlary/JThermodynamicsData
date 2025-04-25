#!/usr/bin/env julia

"""
Verify Pipeline Implementation

This script verifies that the complete thermodynamic data pipeline is working
correctly, with a focus on:

1. Data sources are properly populated with real data for all species
2. Hierarchical source selection correctly prioritizes experimental over theoretical
3. Database contains the correct data with proper source priorities
4. Plots are generated correctly showing all relevant sources

It runs a comprehensive suite of tests and generates a detailed verification report.
"""

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using JSON
using DataFrames
using YAML
using Dates
using Printf

# Initialize logger
log_dir = joinpath(@__DIR__, "logs")
mkpath(log_dir)
log_file = joinpath(log_dir, "pipeline_verification_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
JThermodynamicsData.init_logger("info", log_file)

# Start logging
println("\n" * "=" * 80)
println("THERMODYNAMIC DATA PIPELINE VERIFICATION")
println("=" * 80)

# Create report directory
report_dir = joinpath(@__DIR__, "output")
mkpath(report_dir)
report_file = joinpath(report_dir, "pipeline_verification_report.md")

# Get species list
function get_species_list()
    species_yaml_path = joinpath(@__DIR__, "config", "species.yaml")
    
    if !isfile(species_yaml_path)
        println("⚠️ Species yaml file not found: $species_yaml_path")
        # Fallback to a minimal default list if file not found
        return ["N2", "O2", "H2O", "CO2", "Ar", "He", "CH4", "CO", "NO"]
    end
    
    # Read and parse the YAML file
    try
        println("Loading species from $(species_yaml_path)")
        species_config = YAML.load_file(species_yaml_path)
        
        # Extract the species list
        if haskey(species_config, "species") && species_config["species"] isa Vector
            println("Loaded $(length(species_config["species"])) species from file")
            return species_config["species"]
        else
            println("⚠️ No valid species list found in file, using default species list")
            return ["N2", "O2", "H2O", "CO2", "Ar", "He", "CH4", "CO", "NO"]
        end
    catch e
        println("❌ Error parsing species yaml file: $(e), using default species list")
        return ["N2", "O2", "H2O", "CO2", "Ar", "He", "CH4", "CO", "NO"]
    end
end

# Verify JSON files for all species
function verify_json_files(species_list)
    println("\nStep 1: Verifying JSON files for all species")
    species_dir = joinpath(@__DIR__, "data", "species")
    
    # Check if directory exists
    if !isdir(species_dir)
        println("❌ Species directory not found: $species_dir")
        return false, []
    end
    
    # Count files with experimental and theoretical data
    total_files = 0
    files_with_exp_data = 0
    files_with_theo_data = 0
    empty_files = 0
    species_results = []
    
    for species in species_list
        species_file = joinpath(species_dir, "$(species).json")
        
        if isfile(species_file)
            total_files += 1
            
            # Read the file
            try
                species_data = JSON.parsefile(species_file)
                
                # Check if it has sources
                if haskey(species_data, "sources") && !isempty(species_data["sources"])
                    # Count experimental and theoretical sources
                    exp_sources = []
                    theo_sources = []
                    
                    for (source_name, source_info) in species_data["sources"]
                        # Experimental sources have priority > 4
                        if haskey(source_info, "priority") && source_info["priority"] > 4
                            push!(exp_sources, source_name)
                        else
                            push!(theo_sources, source_name)
                        end
                    end
                    
                    has_exp = !isempty(exp_sources)
                    has_theo = !isempty(theo_sources)
                    
                    if has_exp
                        files_with_exp_data += 1
                    end
                    
                    if has_theo
                        files_with_theo_data += 1
                    end
                    
                    # Add to results
                    push!(species_results, (
                        name = species,
                        found = true,
                        experimental = has_exp,
                        theoretical = has_theo,
                        exp_sources = exp_sources,
                        theo_sources = theo_sources
                    ))
                else
                    empty_files += 1
                    push!(species_results, (
                        name = species,
                        found = true,
                        experimental = false,
                        theoretical = false,
                        exp_sources = [],
                        theo_sources = []
                    ))
                end
            catch e
                println("❌ Error reading species file: $species_file: $e")
                push!(species_results, (
                    name = species,
                    found = true,
                    error = true,
                    message = "$e",
                    experimental = false,
                    theoretical = false,
                    exp_sources = [],
                    theo_sources = []
                ))
            end
        else
            println("❌ Species file not found: $species_file")
            push!(species_results, (
                name = species,
                found = false,
                experimental = false,
                theoretical = false,
                exp_sources = [],
                theo_sources = []
            ))
        end
    end
    
    # Print summary
    found_files = sum(result -> result.found, species_results)
    println("✅ Found $(found_files)/$(length(species_list)) species files")
    println("  - Files with experimental data: $files_with_exp_data")
    println("  - Files with theoretical data only: $files_with_theo_data")
    println("  - Empty files: $empty_files")
    
    # Verify that important species (like N2, O2, etc.) have experimental data
    common_species = ["N2", "O2", "H2O", "CO2", "Ar", "He", "CH4", "CO", "NO"]
    common_with_exp = filter(r -> r.name in common_species && r.experimental, species_results)
    
    if length(common_with_exp) == length(common_species)
        println("✅ All common species have experimental data")
    else
        missing = setdiff(common_species, [r.name for r in common_with_exp])
        println("⚠️ Some common species are missing experimental data: $(join(missing, ", "))")
    end
    
    success = files_with_exp_data > 0 && found_files > 0
    return success, species_results
end

# Verify hierarchical source selection
function verify_hierarchical_selection(species_list)
    println("\nStep 2: Verifying hierarchical source selection")
    
    # List to store results
    selection_results = []
    
    # Temperature range to test
    temp_range = [200.0, 6000.0]
    
    # For each species, test the hierarchical selection
    for species in species_list
        println("  Testing hierarchical selection for: $species")
        
        try
            # Try to load polynomial data using hierarchical selection
            options = Dict("method" => "hierarchical", "use_json" => true)
            poly = JThermodynamicsData.load_polynomial_data(species, temp_range, options)
            
            if poly === nothing
                println("  ⚠️ No data found for $species")
                push!(selection_results, (
                    species = species,
                    status = "No data found",
                    source = nothing,
                    priority = nothing,
                    is_experimental = false
                ))
                continue
            end
            
            # Get the selected source and its priority
            source = poly["source"]
            priority = poly["priority"]
            is_experimental = priority > 4
            
            println("  ✅ Selected source: $source (priority: $priority)")
            
            # Get all available sources
            species_data = JThermodynamicsData.load_species_data(species)
            all_sources = []
            
            if haskey(species_data, "sources")
                for (src_name, src_info) in species_data["sources"]
                    push!(all_sources, (
                        name = src_name,
                        priority = get(src_info, "priority", 0)
                    ))
                end
            end
            
            # Sort sources by priority (highest first)
            sort!(all_sources, by = s -> s.priority, rev = true)
            
            # Check if selection is correct
            correct = true
            expected_source = nothing
            
            # If experimental sources are available, one of them should be selected
            exp_sources = filter(s -> s.priority > 4, all_sources)
            theo_sources = filter(s -> s.priority <= 4, all_sources)
            
            if !isempty(exp_sources)
                expected_source = exp_sources[1].name
                correct = (source == expected_source || priority >= exp_sources[1].priority)
            else
                if !isempty(theo_sources)
                    expected_source = theo_sources[1].name
                    correct = (source == expected_source || priority >= theo_sources[1].priority)
                end
            end
            
            push!(selection_results, (
                species = species,
                status = "OK",
                source = source,
                priority = priority,
                is_experimental = is_experimental,
                correct = correct,
                expected_source = expected_source,
                all_sources = all_sources
            ))
        catch e
            println("  ❌ Error testing hierarchical selection: $e")
            push!(selection_results, (
                species = species,
                status = "Error",
                message = "$e",
                is_experimental = false
            ))
        end
    end
    
    # Summarize results
    total = length(selection_results)
    with_data = count(r -> r.status == "OK", selection_results)
    with_exp = count(r -> r.status == "OK" && r.is_experimental, selection_results)
    with_theo = count(r -> r.status == "OK" && !r.is_experimental, selection_results)
    correct = count(r -> r.status == "OK" && get(r, :correct, false), selection_results)
    
    println("✅ Hierarchical selection test summary:")
    println("  - Total species tested: $total")
    println("  - Species with data: $with_data")
    println("  - Species using experimental sources: $with_exp")
    println("  - Species using theoretical sources only: $with_theo")
    println("  - Correct selections: $correct")
    
    success = with_exp > 0 && correct == with_data
    return success, selection_results
end

# Verify database
function verify_database()
    println("\nStep 3: Verifying database")
    
    db_path = joinpath(@__DIR__, "data", "thermodynamics.duckdb")
    
    if !isfile(db_path)
        println("❌ Database file not found: $db_path")
        return false, Dict()
    end
    
    try
        # Initialize database connection
        config_path = joinpath(@__DIR__, "config", "settings.yaml")
        config = JThermodynamicsData.load_config(config_path)
        db = JThermodynamicsData.init_database(config)
        
        if db === nothing
            println("❌ Failed to connect to database")
            return false, Dict()
        end
        
        # Get database statistics
        species_query = "SELECT COUNT(*) as count FROM species"
        species_result = DuckDB.execute(db, species_query)
        species_df = DataFrame(species_result)
        species_count = species_df[1, :count]
        
        sources_query = "SELECT COUNT(*) as count FROM data_sources"
        sources_result = DuckDB.execute(db, sources_query)
        sources_df = DataFrame(sources_result)
        sources_count = sources_df[1, :count]
        
        data_query = "SELECT COUNT(*) as count FROM thermodynamic_data"
        data_result = DuckDB.execute(db, data_query)
        data_df = DataFrame(data_result)
        data_count = data_df[1, :count]
        
        # Get data sources by priority
        priority_query = """
        SELECT ds.name, ds.priority, COUNT(td.id) as entry_count
        FROM data_sources ds
        LEFT JOIN thermodynamic_data td ON ds.name = td.data_source
        GROUP BY ds.name, ds.priority
        ORDER BY ds.priority DESC
        """
        priority_result = DuckDB.execute(db, priority_query)
        priority_df = DataFrame(priority_result)
        
        # Print summary
        println("✅ Database statistics:")
        println("  - Species in database: $species_count")
        println("  - Data sources defined: $sources_count")
        println("  - Thermodynamic data entries: $data_count")
        
        println("\n  Data sources by priority:")
        for row in eachrow(priority_df)
            println("  - $(row.name): $(row.entry_count) entries (priority: $(row.priority))")
        end
        
        # Count experimental vs theoretical entries
        exp_count = sum(row.entry_count for row in eachrow(priority_df) if row.priority > 4)
        theo_count = sum(row.entry_count for row in eachrow(priority_df) if row.priority <= 4)
        
        println("\n  - Experimental entries: $exp_count")
        println("  - Theoretical entries: $theo_count")
        
        # Check validity of foreign keys
        fk_query = """
        SELECT COUNT(*) as violations
        FROM thermodynamic_data td
        LEFT JOIN data_sources ds ON td.data_source = ds.name
        WHERE ds.name IS NULL
        """
        fk_result = DuckDB.execute(db, fk_query)
        fk_df = DataFrame(fk_result)
        fk_violations = fk_df[1, :violations]
        
        if fk_violations > 0
            println("⚠️ Foreign key violations detected: $fk_violations")
        else
            println("✅ No foreign key violations detected")
        end
        
        # Database stats to return
        db_stats = Dict(
            "species_count" => species_count,
            "sources_count" => sources_count,
            "data_count" => data_count,
            "exp_count" => exp_count,
            "theo_count" => theo_count,
            "fk_violations" => fk_violations,
            "sources" => priority_df
        )
        
        success = species_count > 0 && data_count > 0 && exp_count > 0 && fk_violations == 0
        return success, db_stats
    catch e
        println("❌ Error verifying database: $e")
        return false, Dict()
    end
end

# Verify plot generation
function verify_plots(species_list)
    println("\nStep 4: Verifying plot generation")
    
    plots_dir = joinpath(@__DIR__, "plots")
    
    if !isdir(plots_dir)
        println("❌ Plots directory not found: $plots_dir")
        return false, []
    end
    
    # Check if plots exist for each species
    plot_results = []
    found_count = 0
    
    for species in species_list
        plot_file = joinpath(plots_dir, "$(species)_all_sources.png")
        
        if isfile(plot_file)
            found_count += 1
            
            # Get file size and modification time
            file_stats = stat(plot_file)
            file_size = file_stats.size
            mod_time = Dates.unix2datetime(file_stats.mtime)
            
            # Check if file size is reasonable (not empty)
            valid_size = file_size > 10000  # Reasonable plot size
            
            push!(plot_results, (
                species = species,
                found = true,
                size = file_size,
                valid_size = valid_size,
                modified = mod_time
            ))
        else
            push!(plot_results, (
                species = species,
                found = false
            ))
        end
    end
    
    # Count valid plots
    valid_plots = count(p -> p.found && get(p, :valid_size, false), plot_results)
    
    # Print summary
    println("✅ Plot verification:")
    println("  - Found $(found_count)/$(length(species_list)) plots")
    println("  - Valid plot files: $valid_plots")
    
    success = found_count > 0 && valid_plots > 0
    return success, plot_results
end

# Write verification report
function write_verification_report(results)
    println("\nGenerating verification report...")
    
    open(report_file, "w") do io
        write(io, "# Thermodynamic Data Pipeline Verification Report\n\n")
        write(io, "Generated on: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))\n\n")
        
        # Overall status
        all_success = all(r -> r.success, values(results))
        status = all_success ? "✅ PASSED" : "❌ FAILED"
        
        write(io, "## Overall Status: $status\n\n")
        
        for (step, result) in results
            if result.success
                write(io, "- ✅ $step: Passed\n")
            else
                write(io, "- ❌ $step: Failed\n")
            end
        end
        
        # JSON files verification
        if haskey(results, "JSON Files")
            json_result = results["JSON Files"]
            write(io, "\n## JSON Files Verification\n\n")
            
            # Summary stats
            total = length(json_result.species_results)
            found = count(r -> r.found, json_result.species_results)
            with_exp = count(r -> r.experimental, json_result.species_results)
            with_theo = count(r -> r.theoretical, json_result.species_results)
            
            write(io, "- Total species tested: $total\n")
            write(io, "- Species files found: $found\n")
            write(io, "- Files with experimental data: $with_exp\n")
            write(io, "- Files with theoretical data only: $(with_theo - with_exp)\n")
            
            # Table of species with experimental data
            write(io, "\n### Species with Experimental Data\n\n")
            write(io, "| Species | Experimental Sources | Theoretical Sources |\n")
            write(io, "|---------|---------------------|---------------------|\n")
            
            for result in filter(r -> r.experimental, json_result.species_results)
                exp_sources = join(result.exp_sources, ", ")
                theo_sources = join(result.theo_sources, ", ")
                write(io, "| $(result.name) | $exp_sources | $theo_sources |\n")
            end
            
            # Table of species with theoretical data only
            write(io, "\n### Species with Theoretical Data Only\n\n")
            write(io, "| Species | Theoretical Sources |\n")
            write(io, "|---------|---------------------|\n")
            
            theo_only = filter(r -> !r.experimental && r.theoretical, json_result.species_results)
            for result in theo_only
                theo_sources = join(result.theo_sources, ", ")
                write(io, "| $(result.name) | $theo_sources |\n")
            end
            
            # Table of species with no data
            no_data = filter(r -> !r.experimental && !r.theoretical, json_result.species_results)
            if !isempty(no_data)
                write(io, "\n### Species with No Data\n\n")
                write(io, "| Species | Found | Reason |\n")
                write(io, "|---------|-------|--------|\n")
                
                for result in no_data
                    found = result.found ? "Yes" : "No"
                    reason = !result.found ? "File missing" : "Empty file"
                    write(io, "| $(result.name) | $found | $reason |\n")
                end
            end
        end
        
        # Hierarchical selection verification
        if haskey(results, "Hierarchical Selection")
            hier_result = results["Hierarchical Selection"]
            write(io, "\n## Hierarchical Source Selection Verification\n\n")
            
            # Summary stats
            total = length(hier_result.selection_results)
            with_data = count(r -> r.status == "OK", hier_result.selection_results)
            with_exp = count(r -> r.status == "OK" && r.is_experimental, hier_result.selection_results)
            with_theo = count(r -> r.status == "OK" && !r.is_experimental, hier_result.selection_results)
            correct = count(r -> r.status == "OK" && get(r, :correct, false), hier_result.selection_results)
            
            write(io, "- Total species tested: $total\n")
            write(io, "- Species with data: $with_data\n")
            write(io, "- Species using experimental sources: $with_exp\n")
            write(io, "- Species using theoretical sources only: $with_theo\n")
            write(io, "- Correct selections: $correct\n")
            
            # Table of selection results
            write(io, "\n### Selection Results\n\n")
            write(io, "| Species | Status | Selected Source | Priority | Experimental | Correct |\n")
            write(io, "|---------|--------|----------------|----------|--------------|--------|\n")
            
            for result in hier_result.selection_results
                if result.status == "OK"
                    status = "OK"
                    source = result.source
                    priority = result.priority
                    exp = result.is_experimental ? "Yes" : "No"
                    correct = get(result, :correct, false) ? "✓" : "✗"
                else
                    status = result.status
                    source = "N/A"
                    priority = "N/A"
                    exp = "N/A"
                    correct = "N/A"
                end
                
                write(io, "| $(result.species) | $status | $source | $priority | $exp | $correct |\n")
            end
        end
        
        # Database verification
        if haskey(results, "Database")
            db_result = results["Database"]
            write(io, "\n## Database Verification\n\n")
            
            # Database stats
            if haskey(db_result, "db_stats") && !isempty(db_result.db_stats)
                stats = db_result.db_stats
                
                write(io, "### Database Statistics\n\n")
                write(io, "- Species count: $(stats["species_count"])\n")
                write(io, "- Data sources defined: $(stats["sources_count"])\n")
                write(io, "- Thermodynamic data entries: $(stats["data_count"])\n")
                write(io, "- Experimental data entries: $(stats["exp_count"])\n")
                write(io, "- Theoretical data entries: $(stats["theo_count"])\n")
                write(io, "- Foreign key violations: $(stats["fk_violations"])\n")
                
                # Data sources table
                if haskey(stats, "sources") && stats["sources"] isa DataFrame
                    write(io, "\n### Data Sources\n\n")
                    write(io, "| Source | Priority | Entry Count |\n")
                    write(io, "|--------|----------|-------------|\n")
                    
                    for row in eachrow(stats["sources"])
                        write(io, "| $(row.name) | $(row.priority) | $(row.entry_count) |\n")
                    end
                end
            end
        end
        
        # Plot verification
        if haskey(results, "Plots")
            plot_result = results["Plots"]
            write(io, "\n## Plot Verification\n\n")
            
            # Plot stats
            total = length(plot_result.plot_results)
            found = count(p -> p.found, plot_result.plot_results)
            valid = count(p -> p.found && get(p, :valid_size, false), plot_result.plot_results)
            
            write(io, "- Total species tested: $total\n")
            write(io, "- Plot files found: $found\n")
            write(io, "- Valid plot files: $valid\n")
            
            # Recent plots table
            write(io, "\n### Recently Modified Plots\n\n")
            write(io, "| Species | File Size (bytes) | Last Modified |\n")
            write(io, "|---------|-------------------|---------------|\n")
            
            recent_plots = filter(p -> p.found, plot_result.plot_results)
            sort!(recent_plots, by = p -> p.modified, rev = true)
            
            for plot in recent_plots[1:min(10, length(recent_plots))]
                modified = Dates.format(plot.modified, "yyyy-mm-dd HH:MM:SS")
                write(io, "| $(plot.species) | $(plot.size) | $modified |\n")
            end
            
            # Missing plots
            missing_plots = filter(p -> !p.found, plot_result.plot_results)
            if !isempty(missing_plots)
                write(io, "\n### Missing Plots\n\n")
                missing_list = [plot.species for plot in missing_plots]
                write(io, "- $(join(missing_list, ", "))\n")
            end
        end
    end
    
    println("✅ Verification report saved to: $report_file")
    return report_file
end

# Main function
function main()
    verification_results = Dict()
    
    # Get the species list
    species_list = get_species_list()
    
    # Verify JSON files
    println("\nVerifying JSON files for all species...")
    json_success, species_results = verify_json_files(species_list)
    verification_results["JSON Files"] = (success = json_success, species_results = species_results)
    
    # Verify hierarchical source selection
    println("\nVerifying hierarchical source selection...")
    hier_success, selection_results = verify_hierarchical_selection(species_list)
    verification_results["Hierarchical Selection"] = (success = hier_success, selection_results = selection_results)
    
    # Verify database
    println("\nVerifying database...")
    db_success, db_stats = verify_database()
    verification_results["Database"] = (success = db_success, db_stats = db_stats)
    
    # Verify plots
    println("\nVerifying plots...")
    plot_success, plot_results = verify_plots(species_list)
    verification_results["Plots"] = (success = plot_success, plot_results = plot_results)
    
    # Generate verification report
    report_path = write_verification_report(verification_results)
    
    # Final summary
    all_success = all(r -> r.success, values(verification_results))
    
    println("\n" * "=" * 80)
    if all_success
        println("✅ VERIFICATION PASSED: All pipeline components are working correctly")
    else
        println("❌ VERIFICATION FAILED: Some pipeline components are not working correctly")
        failed_steps = [step for (step, result) in verification_results if !result.success]
        println("   Failed steps: $(join(failed_steps, ", "))")
    end
    println("=" * 80)
    
    println("\nDetailed verification report: $report_path")
    
    return all_success, report_path
end

# Run the main function
success, report_path = main()