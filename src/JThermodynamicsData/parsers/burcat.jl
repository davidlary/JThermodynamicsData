"""
Parser for the Burcat Thermodynamic Database.
Burcat's database is a comprehensive collection of NASA polynomials.
"""

"""
    parse_burcat_file(file_path::String)

Parse a Burcat thermodynamic database file and return a dictionary of species data.
"""
function parse_burcat_file(file_path::String)
    if !isfile(file_path)
        error("Burcat file not found: $file_path")
    end
    
    lines = readlines(file_path)
    return parse_burcat_data(lines)
end

"""
    parse_burcat_data(lines::Vector{String})

Parse Burcat thermodynamic database from a vector of lines.
"""
function parse_burcat_data(lines::Vector{String})
    # Check if the file appears to be HTML content
    if length(lines) > 0 && (
        startswith(lines[1], "<!DOCTYPE") || 
        startswith(lines[1], "<html") || 
        any(l -> occursin("<head>", l) || occursin("<body>", l), lines[1:min(10, length(lines))])
    )
        @error "The file appears to be HTML content rather than a Burcat data file"
        return Dict{String, Dict}()
    end
    
    species_data = Dict{String, Dict}()
    num_parsed = 0
    num_errors = 0
    
    i = 1
    while i <= length(lines)
        line = lines[i]
        
        # Skip empty lines, comments, or HTML content
        if isempty(strip(line)) || startswith(strip(line), "!") || 
           occursin(r"<[a-zA-Z/]", line) || occursin("http://", line) || occursin("https://", line)
            i += 1
            continue
        end
        
        # Burcat format uses a special record structure
        if length(line) >= 80 && !startswith(strip(line), "!")
            # Try to parse a record (4 lines)
            if i + 3 <= length(lines)
                record_lines = lines[i:i+3]
                
                # Check if any of the record lines contain HTML content
                if any(l -> occursin(r"<[a-zA-Z/]", l) || 
                            occursin("http://", l) || 
                            occursin("https://", l), record_lines)
                    i += 1
                    continue
                end
                
                try
                    # Determine if NASA 7 or NASA 9 based on format
                    # Burcat uses a modified NASA 7 format
                    data = parse_burcat_record(record_lines)
                    
                    if !isempty(data)
                        species_name = data["species_name"]
                        species_data[species_name] = data
                        num_parsed += 1
                    end
                catch e
                    @warn "Failed to parse Burcat record at line $i: $e"
                    num_errors += 1
                end
                
                i += 4
            else
                break
            end
        else
            i += 1
        end
    end
    
    @info "Burcat parser: Successfully parsed $num_parsed species with $num_errors errors"
    
    return species_data
end

"""
    parse_burcat_record(lines::Vector{String})

Parse a single Burcat database record (4 lines).
"""
function parse_burcat_record(lines::Vector{String})
    if length(lines) < 4
        error("Invalid Burcat record: not enough lines")
    end
    
    # Extract species information from first line
    header = lines[1]
    
    # Burcat format has more metadata in the first line
    species_name = ""
    formula = ""
    cas_number = ""
    phase = ""
    
    if length(header) >= 18
        species_name = strip(header[1:18])
    end
    
    if length(header) >= 45
        formula = strip(header[19:45])
    end
    
    if length(header) >= 58
        phase = strip(header[46:46])
        cas_number = strip(header[47:58])
    end
    
    # Extract additional metadata
    metadata = Dict{String, Any}(
        "phase" => phase,
        "cas" => cas_number
    )
    
    if length(header) >= 66
        metadata["source"] = strip(header[59:66])
    end
    
    # Parse temperature ranges
    temp_ranges = []
    
    if length(header) >= 105
        try
            t_low = parse(Float64, strip(header[67:75]))
            t_high = parse(Float64, strip(header[76:84]))
            
            # Some Burcat entries use multiple ranges
            if length(header) >= 95 && !isempty(strip(header[85:95]))
                t_mid = parse(Float64, strip(header[85:95]))
                push!(temp_ranges, (t_low, t_mid))
                push!(temp_ranges, (t_mid, t_high))
            else
                push!(temp_ranges, (t_low, t_high))
            end
        catch e
            @warn "Failed to parse temperature ranges: $e"
            # Default range if parsing fails
            push!(temp_ranges, (200.0, 6000.0))
        end
    else
        # Default range if not specified
        push!(temp_ranges, (200.0, 6000.0))
    end
    
    # Molecular weight (if provided)
    mw = 0.0
    if length(header) >= 126 && !isempty(strip(header[106:126]))
        try
            mw = parse(Float64, strip(header[106:126]))
            metadata["molecular_weight"] = mw
        catch e
            @warn "Failed to parse molecular weight: $e"
        end
    end
    
    # Parse coefficients
    # Burcat format has 7 coefficients per temperature range
    n_ranges = length(temp_ranges)
    coefficients = []
    
    # Lines 2-4 contain coefficients
    if n_ranges == 1
        # Single range
        coefs = zeros(7)
        
        if length(lines[2]) >= 75
            for i in 0:4
                start_idx = 1 + i * 15
                end_idx = start_idx + 14
                coefs[i+1] = parse(Float64, strip(lines[2][start_idx:end_idx]))
            end
        end
        
        if length(lines[3]) >= 30
            for i in 0:1
                start_idx = 1 + i * 15
                end_idx = start_idx + 14
                coefs[i+6] = parse(Float64, strip(lines[3][start_idx:end_idx]))
            end
        end
        
        push!(coefficients, coefs)
    else
        # Two ranges - high temperature first in Burcat format
        high_coefs = zeros(7)
        low_coefs = zeros(7)
        
        # High temperature coefficients
        if length(lines[2]) >= 75
            for i in 0:4
                start_idx = 1 + i * 15
                end_idx = start_idx + 14
                high_coefs[i+1] = parse(Float64, strip(lines[2][start_idx:end_idx]))
            end
        end
        
        if length(lines[3]) >= 30
            for i in 0:1
                start_idx = 1 + i * 15
                end_idx = start_idx + 14
                high_coefs[i+6] = parse(Float64, strip(lines[3][start_idx:end_idx]))
            end
        end
        
        # Low temperature coefficients
        if length(lines[3]) >= 75
            for i in 0:2
                start_idx = 31 + i * 15
                end_idx = start_idx + 14
                if end_idx <= length(lines[3])
                    low_coefs[i+1] = parse(Float64, strip(lines[3][start_idx:end_idx]))
                end
            end
        end
        
        if length(lines[4]) >= 75
            for i in 0:3
                start_idx = 1 + i * 15
                end_idx = start_idx + 14
                if i < 2
                    low_coefs[i+4] = parse(Float64, strip(lines[4][start_idx:end_idx]))
                else
                    low_coefs[i+3] = parse(Float64, strip(lines[4][start_idx:end_idx]))
                end
            end
        end
        
        # In Burcat format, the low-temperature range comes first in the temperature order
        push!(coefficients, low_coefs)
        push!(coefficients, high_coefs)
    end
    
    # Calculate molecular weight if not provided
    if mw == 0.0 && !isempty(formula)
        mw = calculate_molecular_weight(formula)
        metadata["molecular_weight"] = mw
    end
    
    return Dict(
        "species_name" => species_name,
        "formula" => formula,
        "temperature_ranges" => temp_ranges,
        "coefficients" => coefficients,
        "polynomial_type" => "nasa7",
        "metadata" => metadata
    )
end

"""
    burcat_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)

Read a Burcat database file and store the data in the database.
"""
function burcat_to_database(conn::DuckDB.DB, file_path::String, source_name::String, reliability_score::Float64=0.0)
    species_data = parse_burcat_file(file_path)
    
    # Transaction for bulk inserts
    transaction(function()
        for (species_name, data) in species_data
            # Get molecular weight from metadata or calculate from formula
            mw = get(data["metadata"], "molecular_weight", 0.0)
            if mw == 0.0 && !isempty(data["formula"])
                mw = calculate_molecular_weight(data["formula"])
            end
            
            # Get CAS number from metadata
            cas = get(data["metadata"], "cas", "")
            
            # Get or create species
            species_id = add_species(conn, species_name, data["formula"], cas, mw, data["metadata"])
            
            # Add thermodynamic data
            temp_ranges = data["temperature_ranges"]
            
            if isempty(temp_ranges)
                @warn "No temperature ranges found for species: $species_name"
                continue
            end
            
            # Temperature range for database (min to max across all ranges)
            temp_min = minimum([range[1] for range in temp_ranges])
            temp_max = maximum([range[2] for range in temp_ranges])
            
            # Create data dictionary
            thermo_data = Dict(
                "coefficients" => data["coefficients"],
                "temperature_ranges" => temp_ranges
            )
            
            # Add to database
            add_thermodynamic_data(
                conn, 
                species_id, 
                source_name, 
                "nasa7", 
                temp_min, 
                temp_max, 
                thermo_data,
                Dict(),  # No uncertainty data
                reliability_score
            )
        end
    end, conn)
    
    return length(species_data)
end

"""
    download_burcat_database(cache_dir::String)

Download the latest Burcat thermodynamic database and save to the cache directory.
"""
function download_burcat_database(cache_dir::String)
    # Primary and backup URLs
    urls = [
        "https://burcat.technion.ac.il/dir/BURCAT.THR", 
        "http://garfield.chem.elte.hu/Burcat/BURCAT.THR"
    ]
    
    output_path = joinpath(cache_dir, "burcat", "BURCAT.THR")
    
    # Create directory if it doesn't exist
    mkpath(dirname(output_path))
    
    for (i, url) in enumerate(urls)
        try
            @info "Trying to download Burcat database from URL $(i) of $(length(urls))"
            
            response = HTTP.get(url, status_exception=false)
            
            if response.status == 200
                content_type = ""
                for header in response.headers
                    if lowercase(header[1]) == "content-type"
                        content_type = lowercase(header[2])
                        break
                    end
                end
                
                # Check if the response is HTML instead of text data
                if contains(content_type, "text/html") || 
                   (length(response.body) > 10 && 
                    (startswith(String(response.body)[1:10], "<!DOCTYPE") || 
                     startswith(String(response.body)[1:10], "<html") ||
                     occursin("<body", String(response.body)[1:min(1000, length(String(response.body)))]))
                    )
                    @warn "URL $(url) returned HTML content instead of data file, trying next URL"
                    continue
                end
                
                open(output_path, "w") do io
                    write(io, response.body)
                end
                
                # Verify the file is a valid Burcat file
                lines = readlines(output_path)
                if !isempty(lines) && !startswith(lines[1], "<!DOCTYPE") && !startswith(lines[1], "<html")
                    valid_lines = 0
                    for line in lines[1:min(10, length(lines))]
                        if length(line) >= 80 && !occursin("<", line) && occursin(r"\d", line)
                            valid_lines += 1
                        end
                    end
                    
                    if valid_lines >= 3
                        @info "Successfully downloaded Burcat database to $(output_path)"
                        return output_path
                    else
                        @warn "Downloaded file doesn't appear to be a valid Burcat database, trying next URL"
                    end
                else
                    @warn "Downloaded file appears to be HTML content, trying next URL"
                end
            else
                @warn "Failed to download Burcat database from $(url): HTTP status $(response.status)"
            end
        catch e
            @warn "Error downloading Burcat database from $(url): $e"
        end
    end
    
    # If all downloads fail, try to use sample data
    sample_path = joinpath(dirname(dirname(dirname(dirname(cache_dir)))), "data", "test_data", "sample_burcat.dat") 
    if isfile(sample_path)
        @info "Using sample Burcat data from $(sample_path)"
        cp(sample_path, output_path, force=true)
        return output_path
    end
    
    error("Failed to download Burcat database from all URLs and no sample data available")
end