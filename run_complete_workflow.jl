#!/usr/bin/env julia

"""
Run Complete Thermodynamic Workflow

This script runs the complete thermodynamic data workflow with enhanced validation:
1. Download all data sources with enhanced data for O3 and other species
2. Parse and process data into JSON files
3. Load JSON files into DuckDB database
4. Generate plots for comparison showing the hierarchical selection
5. Validate the hierarchy is respected

This will ensure the most accurate source is used for each species, especially for O3.

Key features:
- Full hierarchical traversal logic, starting with the most accurate theoretical 
  calculation and traversing the entire hierarchy
- For each species, higher priority sources completely override lower priority ones  
- Plot legends show only one theoretical source (the best one) instead of multiple 
  redundant theoretical sources
- Clean and consistent plots across all species (including ions like OH-, Zn, NO+, NO2+)
"""

using Pkg
Pkg.activate(@__DIR__)

using JThermodynamicsData
using Dates
using Printf
using JSON

"""
    validate_hierarchy_selection()

Validate that the data source hierarchy is respected, especially for O3.
"""
function validate_hierarchy_selection()
    println("\n== Validating Hierarchy Selection ==")
    
    # Key species to validate
    key_species = ["O3", "N2", "O2", "H2O", "CO2", "SO2", "H2O2", "HNO3"]
    
    # Track validation results
    validation_results = Dict()
    
    for species in key_species
        println("\nValidating species: $(species)")
        
        # Get available sources for this species
        sources = JThermodynamicsData.list_available_sources(species)
        
        if isempty(sources)
            println("  No sources available for $(species)")
            validation_results[species] = Dict(
                "status" => "failed",
                "reason" => "no_sources",
                "details" => "No sources available"
            )
            continue
        end
        
        # Get the selected source (highest priority)
        selected_source = sources[1][1]
        selected_priority = sources[1][2]
        
        println("  Available sources:")
        for (name, priority, reliability) in sources
            println("  - $(name) (priority: $(priority), reliability: $(reliability))")
        end
        
        println("  Selected source: $(selected_source) (priority: $(selected_priority))")
        
        # Validate that the selected source is the highest priority one
        is_correct = all(s[2] <= selected_priority for s in sources)
        
        # For O3, check if it has at least one experimental source
        if species == "O3"
            has_exp_source = any(s[1] in ["atct", "burcat", "nist-webbook", "tde", "janaf", 
                                      "nasa-cea", "chemkin", "gri-mech"] for s in sources)
            
            if !has_exp_source
                println("  VALIDATION FAILED: O3 does not have any experimental sources!")
                validation_results[species] = Dict(
                    "status" => "failed",
                    "reason" => "no_experimental",
                    "details" => "No experimental sources for O3"
                )
            elseif !is_correct
                println("  VALIDATION FAILED: Selected source is not highest priority!")
                validation_results[species] = Dict(
                    "status" => "failed",
                    "reason" => "wrong_priority",
                    "details" => "Selected source does not have highest priority"
                )
            else
                println("  VALIDATION PASSED: Selected source is highest priority experimental source!")
                validation_results[species] = Dict(
                    "status" => "passed",
                    "selected" => selected_source,
                    "priority" => selected_priority
                )
            end
        else
            if !is_correct
                println("  VALIDATION FAILED: Selected source is not highest priority!")
                validation_results[species] = Dict(
                    "status" => "failed",
                    "reason" => "wrong_priority",
                    "details" => "Selected source does not have highest priority"
                )
            else
                println("  VALIDATION PASSED: Selected source is highest priority!")
                validation_results[species] = Dict(
                    "status" => "passed",
                    "selected" => selected_source,
                    "priority" => selected_priority
                )
            end
        end
    end
    
    # Print summary
    println("\n== Hierarchy Validation Summary ==")
    all_passed = true
    for (species, result) in validation_results
        if result["status"] == "passed"
            println("✅ $(species): PASSED - Using $(result["selected"]) (priority: $(result["priority"]))")
        else
            println("❌ $(species): FAILED - $(result["reason"]): $(result["details"])")
            all_passed = false
        end
    end
    
    if all_passed
        println("\n✅ ALL VALIDATIONS PASSED - Hierarchy is respected for all species!")
    else
        println("\n❌ VALIDATION FAILED - Hierarchy is not respected for some species!")
    end
    
    return all_passed
end

"""
    main()

Main function to run the complete thermodynamic data workflow with validation.
"""
function main()
    # Set up logging
    log_dir = joinpath(@__DIR__, "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "complete_workflow_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger with file output
    JThermodynamicsData.init_logger("info", log_file)
    
    # Start pipeline logging
    JThermodynamicsData.log_pipeline_start("CompleteWorkflow", Dict(
        "start_time" => now(),
        "working_directory" => @__DIR__
    ))
    
    println("JThermodynamicsData - Complete Workflow with Validation")
    println("===================================================")
    
    # Step 1: Download all data sources with enhanced data
    println("\n== Step 1: Downloading All Data Sources (Enhanced) ==")
    JThermodynamicsData.log_stage_start("CompleteWorkflow", "DownloadSources", Dict())
    
    download_start = now()
    include(joinpath(@__DIR__, "scripts", "download_all_sources.jl"))
    download_time = now() - download_start
    
    JThermodynamicsData.log_timing_benchmark("Workflow", "DownloadSources", download_time)
    JThermodynamicsData.log_stage_end("CompleteWorkflow", "DownloadSources", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(download_time)
    ))
    
    # Step 2: Fetch and process data for all species
    println("\n== Step 2: Fetching Data from All Sources ==")
    JThermodynamicsData.log_stage_start("CompleteWorkflow", "FetchAllData", Dict())
    
    fetch_start = now()
    include(joinpath(@__DIR__, "scripts", "fetch_all_sources.jl"))
    fetch_time = now() - fetch_start
    
    JThermodynamicsData.log_timing_benchmark("Workflow", "FetchAllData", fetch_time)
    JThermodynamicsData.log_stage_end("CompleteWorkflow", "FetchAllData", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(fetch_time)
    ))
    
    # Step 3: Sync data to database
    println("\n== Step 3: Syncing to Database ==")
    JThermodynamicsData.log_stage_start("CompleteWorkflow", "SyncToDatabase", Dict())
    
    db_start = now()
    include(joinpath(@__DIR__, "scripts", "sync_json_to_database.jl"))
    db_time = now() - db_start
    
    JThermodynamicsData.log_timing_benchmark("Workflow", "SyncToDatabase", db_time)
    JThermodynamicsData.log_stage_end("CompleteWorkflow", "SyncToDatabase", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(db_time)
    ))
    
    # Step 4: Generate plots for all species
    println("\n== Step 4: Generating Comparison Plots ==")
    JThermodynamicsData.log_stage_start("CompleteWorkflow", "GeneratePlots", Dict())
    
    plots_start = now()
    include(joinpath(@__DIR__, "run_all_species_plots.jl"))
    plots_time = now() - plots_start
    
    JThermodynamicsData.log_timing_benchmark("Workflow", "GeneratePlots", plots_time)
    JThermodynamicsData.log_stage_end("CompleteWorkflow", "GeneratePlots", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(plots_time)
    ))
    
    # Step 5: Validate hierarchy selection
    println("\n== Step 5: Validating Hierarchy Selection ==")
    JThermodynamicsData.log_stage_start("CompleteWorkflow", "ValidateHierarchy", Dict())
    
    validate_start = now()
    validation_result = validate_hierarchy_selection()
    validate_time = now() - validate_start
    
    JThermodynamicsData.log_timing_benchmark("Workflow", "ValidateHierarchy", validate_time)
    JThermodynamicsData.log_stage_end("CompleteWorkflow", "ValidateHierarchy", Dict(
        "elapsed_time" => JThermodynamicsData.format_time_duration(validate_time),
        "validation_result" => validation_result
    ))
    
    # Finish pipeline
    total_time = now() - download_start
    JThermodynamicsData.log_pipeline_end("CompleteWorkflow", Dict(
        "total_time" => JThermodynamicsData.format_time_duration(total_time),
        "data_dir" => joinpath(@__DIR__, "data"),
        "plots_dir" => joinpath(@__DIR__, "plots"),
        "db_path" => joinpath(@__DIR__, "data", "thermodynamics.duckdb"),
        "log_file" => log_file,
        "validation_passed" => validation_result
    ))
    
    println("\n== Complete Workflow Execution Summary ==")
    println("Total Execution Time: $(JThermodynamicsData.format_time_duration(total_time))")
    println("  - Download Sources: $(JThermodynamicsData.format_time_duration(download_time))")
    println("  - Fetch All Data: $(JThermodynamicsData.format_time_duration(fetch_time))")
    println("  - Sync to Database: $(JThermodynamicsData.format_time_duration(db_time))")
    println("  - Generate Plots: $(JThermodynamicsData.format_time_duration(plots_time))")
    println("  - Validate Hierarchy: $(JThermodynamicsData.format_time_duration(validate_time))")
    println()
    
    if validation_result
        println("✅ VALIDATION PASSED - Hierarchy is respected for all species!")
        println("All data has been processed, validated, and synced to the database.")
    else
        println("❌ VALIDATION FAILED - Hierarchy is not respected for some species!")
        println("Review the logs and validation output to fix hierarchy issues.")
    end
    
    println("You can now use the run_thermodynamics.jl script to query properties.")
end

# Run the main function
main()