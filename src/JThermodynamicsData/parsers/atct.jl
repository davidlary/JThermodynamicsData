"""
Parser for Active Thermochemical Tables (ATcT) data.
ATcT provides highly accurate thermodynamic data based on a thermochemical network approach.
"""

# Reimplement to avoid circular dependencies
function calculate_mw_from_formula(formula::String)
    # ATOMIC_MASSES is defined in core/constants.jl
    elements = Dict{String, Int}()
    
    # Simple regex for formula parsing
    pattern = r"([A-Z][a-z]*)(\d*)"
    matches = eachmatch(pattern, formula)
    
    for m in matches
        element = m.captures[1]
        count = m.captures[2]
        
        if isempty(count)
            count = "1"
        end
        
        elements[element] = get(elements, element, 0) + parse(Int, count)
    end
    
    mw = 0.0
    for (element, count) in elements
        # Use average values for common elements
        if element == "H"
            mw += 1.008 * count
        elseif element == "C"
            mw += 12.011 * count
        elseif element == "N"
            mw += 14.007 * count
        elseif element == "O"
            mw += 15.999 * count
        elseif element == "F"
            mw += 18.998 * count
        elseif element == "Cl"
            mw += 35.453 * count
        elseif element == "Br"
            mw += 79.904 * count
        elseif element == "I"
            mw += 126.904 * count
        elseif element == "S"
            mw += 32.06 * count
        elseif element == "P"
            mw += 30.974 * count
        elseif element == "Si"
            mw += 28.085 * count
        elseif element == "Na"
            mw += 22.99 * count
        elseif element == "K"
            mw += 39.098 * count
        elseif element == "Ca"
            mw += 40.078 * count
        elseif element == "Mg"
            mw += 24.305 * count
        elseif element == "Al"
            mw += 26.982 * count
        elseif element == "Fe"
            mw += 55.845 * count
        elseif element == "Cu"
            mw += 63.546 * count
        elseif element == "Zn"
            mw += 65.38 * count
        elseif element == "Ar"
            mw += 39.948 * count
        elseif element == "He"
            mw += 4.0026 * count
        elseif element == "Ne"
            mw += 20.1797 * count
        elseif element == "Kr"
            mw += 83.798 * count
        elseif element == "Xe"
            mw += 131.293 * count
        else
            @warn "Unknown element in formula: $element"
        end
    end
    
    return mw
end

"""
    parse_atct_file(file_path::String)

Parse an ATcT data file and return a dictionary of thermodynamic data.
"""
function parse_atct_file(file_path::String)
    if !isfile(file_path)
        error("ATcT file not found: $file_path")
    end
    
    lines = readlines(file_path)
    return parse_atct_data(lines)
end

"""
    parse_atct_data(lines::Vector{String})

Parse ATcT data from a vector of lines.
"""
function parse_atct_data(lines::Vector{String})
    species_data = Dict{String, Dict}()
    
    # ATcT files have various formats depending on the version and source
    # This is a generic parser that attempts to handle common formats
    
    i = 1
    current_species = ""
    current_data = nothing
    
    while i <= length(lines)
        line = lines[i]
        stripped = strip(line)
        
        # Skip empty lines and comments
        if isempty(stripped) || startswith(stripped, "#") || startswith(stripped, "//")
            i += 1
            continue
        end
        
        # Check if this is a new species entry
        if occursin("Formula:", line) || (occursin("Species:", line) && !isempty(current_species))
            # Save previous species if we have one
            if !isempty(current_species) && current_data !== nothing
                species_data[current_species] = current_data
            end
            
            # Parse new species
            if occursin("Formula:", line)
                parts = split(line, "Formula:")
                if length(parts) >= 2
                    current_species = strip(parts[2])
                    current_data = Dict(
                        "species_name" => current_species,
                        "formula" => current_species,
                        "metadata" => Dict{String, Any}("source" => "ATcT")
                    )
                end
            elseif occursin("Species:", line)
                parts = split(line, "Species:")
                if length(parts) >= 2
                    current_species = strip(parts[2])
                    current_data = Dict(
                        "species_name" => current_species,
                        "formula" => "",
                        "metadata" => Dict{String, Any}("source" => "ATcT")
                    )
                end
            end
        elseif current_data !== nothing
            # Try to extract formula if not already set
            if isempty(current_data["formula"]) && occursin("Formula:", line)
                parts = split(line, "Formula:")
                if length(parts) >= 2
                    current_data["formula"] = strip(parts[2])
                end
            end
            
            # Extract CAS number
            if occursin("CAS:", line)
                parts = split(line, "CAS:")
                if length(parts) >= 2
                    current_data["metadata"]["cas"] = strip(parts[2])
                end
            end
            
            # Extract standard enthalpy of formation
            if occursin("ΔfH°(298.15 K):", line) || occursin("Enthalpy of formation:", line)
                parts = split(line, ":")
                if length(parts) >= 2
                    value_part = strip(parts[2])
                    
                    # Extract value and uncertainty
                    value_match = match(r"([\d\.\-]+)\s*±\s*([\d\.]+)", value_part)
                    if value_match !== nothing
                        value = parse(Float64, value_match.captures[1])
                        uncertainty = parse(Float64, value_match.captures[2])
                        
                        if !haskey(current_data, "properties")
                            current_data["properties"] = Dict{String, Any}()
                        end
                        
                        current_data["properties"]["enthalpy_formation"] = Dict(
                            "value" => value,
                            "uncertainty" => uncertainty,
                            "temperature" => 298.15
                        )
                    else
                        # Try to parse just the value
                        value_match = match(r"([\d\.\-]+)", value_part)
                        if value_match !== nothing
                            value = parse(Float64, value_match.captures[1])
                            
                            if !haskey(current_data, "properties")
                                current_data["properties"] = Dict{String, Any}()
                            end
                            
                            current_data["properties"]["enthalpy_formation"] = Dict(
                                "value" => value,
                                "temperature" => 298.15
                            )
                        end
                    end
                end
            end
            
            # Extract standard entropy
            if occursin("S°(298.15 K):", line) || occursin("Entropy:", line)
                parts = split(line, ":")
                if length(parts) >= 2
                    value_part = strip(parts[2])
                    
                    # Extract value and uncertainty
                    value_match = match(r"([\d\.\-]+)\s*±\s*([\d\.]+)", value_part)
                    if value_match !== nothing
                        value = parse(Float64, value_match.captures[1])
                        uncertainty = parse(Float64, value_match.captures[2])
                        
                        if !haskey(current_data, "properties")
                            current_data["properties"] = Dict{String, Any}()
                        end
                        
                        current_data["properties"]["entropy"] = Dict(
                            "value" => value,
                            "uncertainty" => uncertainty,
                            "temperature" => 298.15
                        )
                    else
                        # Try to parse just the value
                        value_match = match(r"([\d\.\-]+)", value_part)
                        if value_match !== nothing
                            value = parse(Float64, value_match.captures[1])
                            
                            if !haskey(current_data, "properties")
                                current_data["properties"] = Dict{String, Any}()
                            end
                            
                            current_data["properties"]["entropy"] = Dict(
                                "value" => value,
                                "temperature" => 298.15
                            )
                        end
                    end
                end
            end
            
            # Extract heat capacity
            if occursin("Cp(298.15 K):", line) || occursin("Heat capacity:", line)
                parts = split(line, ":")
                if length(parts) >= 2
                    value_part = strip(parts[2])
                    
                    # Extract value and uncertainty
                    value_match = match(r"([\d\.\-]+)\s*±\s*([\d\.]+)", value_part)
                    if value_match !== nothing
                        value = parse(Float64, value_match.captures[1])
                        uncertainty = parse(Float64, value_match.captures[2])
                        
                        if !haskey(current_data, "properties")
                            current_data["properties"] = Dict{String, Any}()
                        end
                        
                        current_data["properties"]["heat_capacity"] = Dict(
                            "value" => value,
                            "uncertainty" => uncertainty,
                            "temperature" => 298.15
                        )
                    else
                        # Try to parse just the value
                        value_match = match(r"([\d\.\-]+)", value_part)
                        if value_match !== nothing
                            value = parse(Float64, value_match.captures[1])
                            
                            if !haskey(current_data, "properties")
                                current_data["properties"] = Dict{String, Any}()
                            end
                            
                            current_data["properties"]["heat_capacity"] = Dict(
                                "value" => value,
                                "temperature" => 298.15
                            )
                        end
                    end
                end
            end
            
            # Extract tabular data if in a table format
            if occursin("T(K)", line) && occursin("Cp", line)
                # This might be the start of a table with temperature-dependent properties
                table_headers = split(stripped)
                
                # Create table data structure
                tabular_data = Dict{String, Vector{Float64}}()
                for header in table_headers
                    tabular_data[header] = Float64[]
                end
                
                # Read table rows
                i += 1
                while i <= length(lines)
                    data_line = strip(lines[i])
                    
                    # Check if we've reached the end of the table
                    if isempty(data_line) || !any(isdigit, data_line)
                        break
                    end
                    
                    # Parse data row
                    data_values = split(data_line)
                    if length(data_values) == length(table_headers)
                        for j in 1:length(table_headers)
                            try
                                value = parse(Float64, data_values[j])
                                push!(tabular_data[table_headers[j]], value)
                            catch e
                                @warn "Failed to parse table value: $(data_values[j])"
                            end
                        end
                    end
                    
                    i += 1
                end
                
                # Set tabular data if we found any
                if !isempty(tabular_data) && haskey(tabular_data, "T(K)") && !isempty(tabular_data["T(K)"])
                    current_data["tabular_data"] = tabular_data
                    
                    # Get temperature range
                    temps = tabular_data["T(K)"]
                    temp_min = minimum(temps)
                    temp_max = maximum(temps)
                    current_data["temperature_range"] = [temp_min, temp_max]
                    
                    # Continue to next line (already incremented i)
                    continue
                end
            end
        end
        
        i += 1
    end
    
    # Save last species if we have one
    if !isempty(current_species) && current_data !== nothing
        species_data[current_species] = current_data
    end
    
    return species_data
end

"""
    atct_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read an ATcT data file and store the data in the database.
"""
function atct_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_atct_file(file_path)
    
    count = 0
    
    # Transaction for bulk inserts
    transaction(function()
        for (species_name, data) in species_data
            # Calculate molecular weight if formula is available
            mw = 0.0
            formula = data["formula"]
            
            if !isempty(formula)
                mw = calculate_mw_from_formula(formula)
                data["metadata"]["molecular_weight"] = mw
            end
            
            # Get CAS number if available
            cas = get(data["metadata"], "cas", "")
            
            # Get or create species
            species_id = add_species(conn, species_name, formula, cas, mw, data["metadata"])
            
            # Add thermodynamic data
            if haskey(data, "tabular_data")
                # Use tabular data if available
                temp_range = data["temperature_range"]
                temp_min = temp_range[1]
                temp_max = temp_range[2]
                
                # Create data dictionary
                thermo_data = Dict(
                    "tabular_data" => data["tabular_data"]
                )
                
                # Create uncertainty dictionary
                uncertainty = Dict{String, Any}()
                
                # Add to database
                add_thermodynamic_data(
                    conn, 
                    species_id, 
                    source_name, 
                    "tabular", 
                    temp_min, 
                    temp_max, 
                    thermo_data,
                    uncertainty,
                    reliability_score
                )
            elseif haskey(data, "properties")
                # Use single-point properties
                props = data["properties"]
                
                # Set temperature range to standard temperature
                temp_min = 298.15
                temp_max = 298.15
                
                # Create data dictionary
                thermo_data = Dict(
                    "properties" => props
                )
                
                # Create uncertainty dictionary
                uncertainty = Dict{String, Any}()
                
                # Extract uncertainty values
                for (prop_name, prop_data) in props
                    if haskey(prop_data, "uncertainty")
                        uncertainty[prop_name] = prop_data["uncertainty"]
                    end
                end
                
                # Add to database
                add_thermodynamic_data(
                    conn, 
                    species_id, 
                    source_name, 
                    "properties", 
                    temp_min, 
                    temp_max, 
                    thermo_data,
                    uncertainty,
                    reliability_score
                )
            else
                # Skip if no usable data
                continue
            end
            
            count += 1
        end
    end, conn)
    
    return count
end

"""
    download_atct_data(cache_dir::String)

Download the latest ATcT data from the Argonne National Laboratory website.
"""
function download_atct_data(cache_dir::String)
    # ATcT data is typically available as special files
    # Here we're using a simplified approach that would need to be updated
    # with the actual URLs for the ATcT data files
    
    url = "https://atct.anl.gov/Thermochemical%20Data/version%201.122/ATcT.nist"
    output_path = joinpath(cache_dir, "atct", "ATcT.nist")
    
    # Create directory if it doesn't exist
    mkpath(dirname(output_path))
    
    try
        response = HTTP.get(url)
        
        if response.status == 200
            open(output_path, "w") do io
                write(io, response.body)
            end
            return output_path
        else
            error("Failed to download ATcT data: HTTP status $(response.status)")
        end
    catch e
        error("Failed to download ATcT data: $e")
    end
end