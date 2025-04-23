#!/usr/bin/env julia

"""
Test script to demonstrate the enhanced logging functionality.
"""

using Pkg
Pkg.activate(dirname(dirname(@__FILE__)))

using JThermodynamicsData
using Dates

function test_logging()
    # Set up a log file
    log_dir = joinpath(dirname(dirname(@__FILE__)), "logs")
    mkpath(log_dir)
    log_file = joinpath(log_dir, "logging_test_$(Dates.format(now(), "yyyymmdd_HHMMSS")).log")
    
    # Initialize logger
    JThermodynamicsData.init_logger("debug", log_file)
    
    # Log pipeline start
    JThermodynamicsData.log_pipeline_start("LoggingTest", Dict(
        "test_id" => rand(1000:9999),
        "timestamp" => now()
    ))
    
    # Simulate a multi-stage pipeline
    for stage_num in 1:3
        stage_name = "Stage$(stage_num)"
        
        # Log stage start
        JThermodynamicsData.log_stage_start("LoggingTest", stage_name, Dict(
            "stage_number" => stage_num,
            "started_at" => now()
        ))
        
        # Simulate some work with timing
        JThermodynamicsData.with_timing(
            () -> begin
                # Simulate work
                sleep(0.5 + rand())
                # Log some operations within the stage
                JThermodynamicsData.log_operation("ProcessData", Dict(
                    "data_points" => rand(100:500),
                    "algorithm" => "test-algorithm-$(stage_num)"
                ))
                sleep(0.2)
            end,
            "Testing",
            "Stage$(stage_num)Work"
        )
        
        # Log stage end
        JThermodynamicsData.log_stage_end("LoggingTest", stage_name, Dict(
            "result" => "success",
            "details" => "Completed test stage $(stage_num)"
        ))
    end
    
    # Test error logging
    try
        error("This is a test error")
    catch e
        JThermodynamicsData.log_error("Encountered test error", Dict(
            "error" => e,
            "context" => "error testing phase"
        ))
    end
    
    # Simulate species processing
    JThermodynamicsData.log_species_processing("TestSpecies", [300.0, 1000.0], 10.0)
    
    # Test parser logging
    JThermodynamicsData.log_parser_operation("TestParser", "JANAF", "/path/to/testfile.dat")
    
    # Finish pipeline
    JThermodynamicsData.log_pipeline_end("LoggingTest", Dict(
        "status" => "complete",
        "stages_completed" => 3
    ))
    
    @info "Logging test complete. Log file written to $log_file"
    
    return log_file
end

# Run the test
log_file = test_logging()
println("Log file created at: $log_file")
println("You can view it with: cat $log_file")