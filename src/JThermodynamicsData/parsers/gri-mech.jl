module JThermodynamicsData

export parse_gri_mech

using Printf

"""
    parse_gri_mech(species_name::String; debug::Bool=false)

Parse thermodynamic data for a species from GRI-Mech 3.0 database.
GRI-Mech provides thermodynamic data for natural gas combustion in CHEMKIN format.

# Arguments
- `species_name::String`: Name of the species to fetch
- `debug::Bool=false`: Enable debug output

# Returns
- `Dict`: Thermodynamic data in standardized format
"""
function parse_gri_mech(species_name::String; debug::Bool=false)
    if debug
        println("Parsing GRI-Mech data for $species_name")
    end
    
    # Initialize result structure
    result = Dict{String,Any}(
        "polynomial_type" => "nasa7",
        "temperature_min" => 0.0,
        "temperature_max" => 0.0,
        "data" => Dict{String,Any}(),
        "uncertainty" => 0.03  # Typical uncertainty for GRI-Mech data
    )
    
    # Path to local GRI-Mech data file
    data_path = joinpath(dirname(dirname(dirname(@__FILE__))), "data", "cache", "gri-mech", "therm.dat")
    
    if !isfile(data_path)
        if debug
            println("GRI-Mech data file not found at $data_path")
        end
        return result
    end
    
    try
        # GRI-Mech uses the same format as CHEMKIN, so we can use the same parsing logic
        
        # Read the GRI-Mech data file
        lines = readlines(data_path)
        
        # Check for metadata lines at the beginning
        line_index = 1
        
        # Skip any header comments (lines starting with '!')
        while line_index <= length(lines) && (isempty(strip(lines[line_index])) || startswith(strip(lines[line_index]), "!"))
            line_index += 1
        end
        
        # Check for THERMO line
        if line_index <= length(lines) && occursin("THERMO", uppercase(lines[line_index]))
            line_index += 1
        end
        
        # Check for temperature ranges line
        temp_ranges = [0.0, 0.0, 0.0]  # Default values
        
        if line_index <= length(lines)
            temp_line = lines[line_index]
            
            # Try to parse temperature ranges
            if length(temp_line) >= 30 && all(isdigit, filter(c -> !isspace(c), temp_line[1:30]))
                temp_ranges = [
                    parse(Float64, strip(temp_line[1:10])),
                    parse(Float64, strip(temp_line[11:20])),
                    parse(Float64, strip(temp_line[21:30]))
                ]
                line_index += 1
            end
        end
        
        # Search for the species
        while line_index <= length(lines)
            # GRI-Mech/CHEMKIN format has 4 lines per species
            # Line 1: Species name, elemental composition, temperature ranges, etc.
            # Line 2: Coefficients a1-a5 for high temperature range
            # Line 3: Coefficients a6-a7 for high temperature range, a1-a3 for low temperature range
            # Line 4: Coefficients a4-a7 for low temperature range
            
            # Check if we have enough lines left
            if line_index + 3 > length(lines)
                break
            end
            
            # Check if the current line contains the species name
            current_line = lines[line_index]
            
            if length(current_line) >= 18
                line_species = strip(current_line[1:18])
                
                # GRI-Mech often uses different naming conventions
                # E.g., CH4 might be labeled as "CH4" or "CH4    "
                if lowercase(line_species) == lowercase(species_name) || 
                   lowercase(strip(line_species)) == lowercase(species_name)
                    # Found the species, extract data
                    if debug
                        println("Found species $species_name at line $line_index")
                    end
                    
                    # Extract species-specific temperature ranges if available
                    species_temp_ranges = temp_ranges
                    
                    if length(current_line) >= 75
                        try
                            species_temp_ranges = [
                                parse(Float64, strip(current_line[46:55])),
                                parse(Float64, strip(current_line[56:65])),
                                parse(Float64, strip(current_line[66:75]))
                            ]
                        catch
                            # Use default temperature ranges
                        end
                    end
                    
                    # Extract coefficients
                    line1 = current_line
                    line2 = lines[line_index + 1]
                    line3 = lines[line_index + 2]
                    line4 = lines[line_index + 3]
                    
                    # Extract high temperature coefficients
                    high_coeffs = extract_gri_mech_coefficients([line2, line3[1:30]], debug)
                    
                    # Extract low temperature coefficients
                    low_coeffs = extract_gri_mech_coefficients([line3[31:end], line4], debug)
                    
                    # Set temperature ranges
                    result["temperature_min"] = species_temp_ranges[1]
                    result["temperature_max"] = species_temp_ranges[3]
                    
                    # Store coefficients
                    result["data"]["low_temp"] = Dict(
                        "range" => [species_temp_ranges[1], species_temp_ranges[2]],
                        "coefficients" => low_coeffs
                    )
                    
                    result["data"]["high_temp"] = Dict(
                        "range" => [species_temp_ranges[2], species_temp_ranges[3]],
                        "coefficients" => high_coeffs
                    )
                    
                    break
                end
            end
            
            # Move to the next species
            line_index += 4
        end
        
        if isempty(result["data"])
            if debug
                println("Species $species_name not found in GRI-Mech data file")
            end
        end
    catch e
        if debug
            println("Error parsing GRI-Mech data: $e")
        end
    end
    
    # Return the parsed data
    return result
end

"""
    extract_gri_mech_coefficients(lines::Vector{String}, debug::Bool)

Extract coefficients from GRI-Mech format lines.

# Arguments
- `lines::Vector{String}`: Lines containing GRI-Mech format coefficients
- `debug::Bool`: Enable debug output

# Returns
- `Vector{Float64}`: Array of coefficients
"""
function extract_gri_mech_coefficients(lines::Vector{String}, debug::Bool)
    coeffs = Float64[]
    
    try
        # Each line contains up to 5 coefficients, 15 characters each
        for line in lines
            # Skip empty lines
            if isempty(strip(line))
                continue
            end
            
            # Parse coefficients
            for i in 1:5
                start_idx = (i - 1) * 15 + 1
                end_idx = min(start_idx + 14, length(line))
                
                if end_idx < start_idx
                    # Line too short
                    continue
                end
                
                coeff_str = line[start_idx:end_idx]
                
                # Skip empty coefficient slots
                if isempty(strip(coeff_str))
                    continue
                end
                
                # Replace D or d with E for scientific notation
                coeff_str = replace(coeff_str, "D" => "E")
                coeff_str = replace(coeff_str, "d" => "E")
                
                # Parse coefficient
                try
                    coeff = parse(Float64, coeff_str)
                    push!(coeffs, coeff)
                catch e
                    if debug
                        println("Error parsing coefficient '$coeff_str': $e")
                    end
                    push!(coeffs, 0.0)
                end
            end
        end
        
        # Ensure we have exactly 7 coefficients
        while length(coeffs) < 7
            push!(coeffs, 0.0)
        end
        
        if length(coeffs) > 7
            coeffs = coeffs[1:7]
        end
    catch e
        if debug
            println("Error extracting GRI-Mech coefficients: $e")
        end
        
        # Return default coefficients on error
        return [3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    end
    
    return coeffs
end

end # module
