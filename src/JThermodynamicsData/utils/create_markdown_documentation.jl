function create_markdown_documentation(species_name::String, result::Dict, all_sources::Vector{Dict}, output_dir::String)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Get formula
    formula = result["formula"]
    temperature = result["temperature"]
    
    # Sort sources by priority
    # First theoretical, then experimental by data source priority
    sorted_sources = sort(all_sources, by = s -> startswith(s["data_source"], "THEORETICAL") ? 0 : 
                          (haskey(s, "priority") ? s["priority"] : 9))
    
    # Determine the final source (highest priority source)
    final_source = sorted_sources[end]["data_source"]
    final_priority = haskey(sorted_sources[end], "priority") ? sorted_sources[end]["priority"] : 0
    
    # Count theoretical vs experimental sources
    theoretical_count = count(s -> startswith(s["data_source"], "THEORETICAL"), sorted_sources)
    experimental_count = length(sorted_sources) - theoretical_count
    
    # Create markdown file
    markdown_file = joinpath(output_dir, "$(species_name)_sources.md")
    open(markdown_file, "w") do io
        write(io, "# Thermodynamic Data Sources for $species_name ($formula)\n\n")
        
        # Summary
        write(io, "## Summary\n\n")
        write(io, "- Temperature: $(temperature) K\n")
        write(io, "- Theoretical methods used: $theoretical_count\n")
        write(io, "- Experimental sources used: $experimental_count\n")
        write(io, "- Final values taken from: $final_source (Priority: $final_priority)\n\n")
        
        # Final properties
        write(io, "## Final Properties\n\n")
        write(io, "| Property | Value | Uncertainty | Units | Source |\n")
        write(io, "|----------|-------|------------|-------|--------|\n")
        
        for prop in ["Cp", "H", "S", "G"]
            if haskey(result["properties"], prop)
                value = result["properties"][prop]["value"]
                uncertainty = result["properties"][prop]["uncertainty"]
                units = result["properties"][prop]["units"]
                write(io, "| $prop | $value | $uncertainty | $units | $final_source |\n")
            end
        end
        
        write(io, "\n")
        
        # Source details
        write(io, "## Data Sources\n\n")
        write(io, "**Note: The final refined values are taken directly from the highest priority source available.**\n")
        write(io, "Uncertainty is calculated from the spread of values across all available sources.\n\n")
        write(io, "Sources listed in order of priority, from theoretical methods to highest priority experimental:\n\n")
        
        # Create table of sources
        write(io, "| Source | Type | Priority | Used for Values | Used for Uncertainty |\n")
        write(io, "|--------|------|----------|----------------|---------------------|\n")
        
        for (i, source) in enumerate(sorted_sources)
            # Determine source type
            source_type = startswith(source["data_source"], "THEORETICAL") ? "Theoretical" : "Experimental"
            
            # Determine priority
            source_priority = haskey(source, "priority") ? source["priority"] : 0
            
            # Check if this source was used for the final values
            used_for_values = (i == length(sorted_sources)) ? "Yes" : "No"
            
            # All sources contribute to uncertainty
            used_for_uncertainty = "Yes"
            
            write(io, "| $(source["data_source"]) | $source_type | $source_priority | $used_for_values | $used_for_uncertainty |\n")
        end
        
        write(io, "\n")
        
        # Add property values from each source
        write(io, "## Property Values by Source\n\n")
        write(io, "This table shows values from each source, illustrating how uncertainty is derived from the spread of values:\n\n")
        
        # Table headers for each property
        write(io, "| Source | ")
        for prop in ["Cp", "H", "S", "G"]
            write(io, "$prop | ")
        end
        write(io, "\n")
        
        # Separator line
        write(io, "|--------|")
        for _ in 1:4
            write(io, "------|")
        end
        write(io, "\n")
        
        # Values from each source
        for source in sorted_sources
            write(io, "| $(source["data_source"]) | ")
            
            for prop in ["Cp", "H", "S", "G"]
                if haskey(source["properties"], prop)
                    value = source["properties"][prop]["value"]
                    write(io, "$(round(value, digits=4)) | ")
                else
                    write(io, "N/A | ")
                end
            end
            write(io, "\n")
        end
        
        write(io, "\n")
        
        # Add metadata
        write(io, "## Metadata\n\n")
        write(io, "- **Generated**: $(Dates.now())\n")
        write(io, "- **JThermodynamicsData Version**: 1.0.0\n")
        
        # Add polynomial information if available
        if haskey(result, "polynomials")
            write(io, "- **NASA-7 Polynomial File**: $(species_name)_nasa7.txt\n")
            write(io, "- **NASA-9 Polynomial File**: $(species_name)_nasa9.txt\n")
        end
    end
    
    return markdown_file
end
EOL < /dev/null