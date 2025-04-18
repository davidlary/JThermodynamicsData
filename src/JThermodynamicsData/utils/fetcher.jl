"""
    JThermodynamicsData.utils.fetcher

This module handles fetching data from various thermodynamic data sources.
- Retrieves data from original sources
- Implements local caching
- Checks for newer versions when cached files exist
"""

using HTTP
using Dates
using JLD2
using JSON
using YAML
using DataFrames
using CSV

"""
    fetch_data_source(source_name::String, options::Dict=Dict())

Fetch data from a specific source with caching support.
Handles downloading, parsing, and caching the data.
"""
function fetch_data_source(source_name::String, options::Dict=Dict())
    # Default options
    cache_dir = get(options, "cache_dir", joinpath("data", "cache", lowercase(source_name)))
    cache_file = get(options, "cache_file", joinpath(cache_dir, "$(lowercase(source_name)).jld2"))
    force_refresh = get(options, "force_refresh", false)
    refresh_interval = get(options, "refresh_interval", 7)  # days
    
    # Create cache directory if needed
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    
    # Check for cache and use it if valid
    if !force_refresh && isfile(cache_file)
        cache_info = check_cache(cache_file, refresh_interval)
        if cache_info["valid"]
            @info "Using cached $(source_name) data from $(cache_info["date"])"
            return load(cache_file, "data")
        end
    end
    
    # No valid cache, fetch from source
    @info "Fetching $(source_name) data from original source..."
    
    # Dispatch to appropriate fetcher based on source name
    data = Dict()
    
    if lowercase(source_name) == "burcat"
        data = fetch_burcat(options)
    elseif lowercase(source_name) == "nasa"
        data = fetch_nasa(options)
    elseif lowercase(source_name) == "janaf"
        data = fetch_janaf(options)
    elseif lowercase(source_name) == "nist"
        data = fetch_nist(options)
    elseif lowercase(source_name) == "atct"
        data = fetch_atct(options)
    elseif lowercase(source_name) == "tde"
        data = fetch_tde(options)
    elseif lowercase(source_name) == "thermoml"
        data = fetch_thermoml(options)
    elseif lowercase(source_name) == "gri-mech"
        data = fetch_gri_mech(options)
    elseif lowercase(source_name) == "cas"
        data = fetch_cas_registry(options)
    else
        error("Unknown data source: $(source_name)")
    end
    
    # Add metadata
    data["source_name"] = source_name
    data["download_date"] = now()
    data["version"] = get(data, "version", "1.0.0")
    
    # Cache the data
    update_cache(cache_file, data)
    
    return data
end

"""
    check_cache(cache_file::String, refresh_interval::Int=7)

Check if a cache file is valid and within the refresh interval.
"""
function check_cache(cache_file::String, refresh_interval::Int=7)
    result = Dict(
        "valid" => false,
        "date" => nothing,
        "reason" => "Cache file does not exist"
    )
    
    if !isfile(cache_file)
        return result
    end
    
    try
        # Load metadata from cache
        file = jldopen(cache_file, "r")
        if !haskey(file, "metadata")
            close(file)
            return Dict("valid" => false, "reason" => "No metadata in cache")
        end
        
        metadata = file["metadata"]
        close(file)
        
        # Check timestamp
        if !haskey(metadata, "date")
            return Dict("valid" => false, "reason" => "No date in metadata")
        end
        
        cache_date = metadata["date"]
        result["date"] = cache_date
        
        # Check if cache is too old
        days_since_update = (now() - cache_date).value / (1000 * 60 * 60 * 24)
        if days_since_update > refresh_interval
            return Dict(
                "valid" => false, 
                "date" => cache_date,
                "reason" => "Cache is $(round(days_since_update, digits=1)) days old (refresh interval: $(refresh_interval) days)"
            )
        end
        
        # Check version if available
        if haskey(metadata, "source_version") && haskey(metadata, "latest_check")
            latest_check = metadata["latest_check"]
            days_since_check = (now() - latest_check).value / (1000 * 60 * 60 * 24)
            
            # Check for newer version once a day at most
            if days_since_check >= 1
                source_name = metadata["source_name"]
                current_version = metadata["source_version"]
                latest_version = check_latest_version(source_name)
                
                if latest_version != nothing && latest_version != current_version
                    return Dict(
                        "valid" => false, 
                        "date" => cache_date,
                        "reason" => "Newer version available: $(latest_version) (current: $(current_version))"
                    )
                end
                
                # Update last check time
                update_metadata(cache_file, Dict("latest_check" => now()))
            end
        end
        
        # All checks passed
        return Dict("valid" => true, "date" => cache_date)
        
    catch e
        return Dict("valid" => false, "reason" => "Error checking cache: $(e)")
    end
end

"""
    update_cache(cache_file::String, data::Dict)

Update the cache file with new data and metadata.
"""
function update_cache(cache_file::String, data::Dict)
    try
        # Create metadata
        metadata = Dict(
            "date" => get(data, "download_date", now()),
            "source_name" => get(data, "source_name", "unknown"),
            "source_version" => get(data, "version", "1.0.0"),
            "latest_check" => now()
        )
        
        # Save data and metadata
        jldopen(cache_file, "w") do file
            file["data"] = data
            file["metadata"] = metadata
        end
        
        @info "Updated cache file: $(cache_file)"
        return true
    catch e
        @error "Failed to update cache: $(e)"
        return false
    end
end

"""
    update_metadata(cache_file::String, updates::Dict)

Update specific metadata fields in the cache file.
"""
function update_metadata(cache_file::String, updates::Dict)
    try
        file = jldopen(cache_file, "r+")
        if !haskey(file, "metadata")
            metadata = Dict()
        else
            metadata = file["metadata"]
        end
        
        # Update fields
        for (key, value) in updates
            metadata[key] = value
        end
        
        file["metadata"] = metadata
        close(file)
        
        return true
    catch e
        @error "Failed to update metadata: $(e)"
        return false
    end
end

"""
    check_latest_version(source_name::String)

Check for the latest version of a data source.
Returns the version string if available, or nothing if can't determine.
"""
function check_latest_version(source_name::String)
    try
        if lowercase(source_name) == "burcat"
            # Check Burcat database version (from last update date on website)
            response = HTTP.get("http://garfield.chem.elte.hu/Burcat/THERM.DAT")
            if response.status == 200
                # Extract date from header or footer
                content = String(response.body)
                # Simple heuristic: use last modified date as version
                return "$(Dates.today())"
            end
        elseif lowercase(source_name) == "nasa"
            # NASA CEA typically doesn't have version in file, use HTTP last-modified
            response = HTTP.head("https://cearun.grc.nasa.gov/cea/thermo.inp")
            if haskey(response.headers, "Last-Modified")
                last_mod = response.headers["Last-Modified"]
                try
                    date = Dates.DateTime(last_mod, "e, d u Y H:M:S GMT")
                    return "$(date)"
                catch
                    return "$(Dates.today())"
                end
            end
        end
        # For other sources, implement similar version checking logic
        
        # Default behavior if no specific check is implemented
        return nothing
    catch e
        @warn "Error checking latest version for $(source_name): $(e)"
        return nothing
    end
end

"""
    fetch_cas_registry(options::Dict=Dict())

Fetch CAS Registry Numbers for chemical species.
"""
function fetch_cas_registry(options::Dict=Dict())
    cache_dir = get(options, "cache_dir", joinpath("data", "cache", "cas"))
    cache_file = get(options, "cache_file", joinpath(cache_dir, "cas_registry.jld2"))
    
    # Create cache directory if needed
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    
    # Use cached data if available and not forced to refresh
    force_refresh = get(options, "force_refresh", false)
    if !force_refresh && isfile(cache_file)
        cache_info = check_cache(cache_file)
        if cache_info["valid"]
            @info "Using cached CAS registry data"
            return load(cache_file, "data")
        end
    end
    
    # Fallback to built-in data if can't download
    cas_data = Dict(
        "species_to_cas" => Dict(
            # Common species
            "N2" => "7727-37-9",
            "O2" => "7782-44-7",
            "H2" => "1333-74-0", 
            "H2O" => "7732-18-5",
            "CO2" => "124-38-9",
            "CO" => "630-08-0",
            "CH4" => "74-82-8",
            "Ar" => "7440-37-1",
            "He" => "7440-59-7",
            "Ne" => "7440-01-9",
            "Xe" => "7440-63-3",
            "NO" => "10102-43-9",
            "NO2" => "10102-44-0",
            "N2O" => "10024-97-2",
            "O3" => "10028-15-6",
            "OH" => "3352-57-6",
            "H2O2" => "7722-84-1",
            
            # Ions
            "H+" => "12408-02-5",
            "O+" => "17778-80-2",
            "N+" => "14701-21-0",
            "C+" => "14067-05-3",
            "He+" => "14802-31-2",
            "Ar+" => "14791-70-5",
            "e-" => "15380-73-3"
        ),
        "cas_to_species" => Dict(),  # Will be populated below
        "source_name" => "CAS",
        "download_date" => now(),
        "version" => "2023-02-01"
    )
    
    # Create reverse mapping
    for (species, cas) in cas_data["species_to_cas"]
        cas_data["cas_to_species"][cas] = species
    end
    
    # Cache the data
    update_cache(cache_file, cas_data)
    
    return cas_data
end

# Specific fetchers for each data source
# These will be implemented with the appropriate logic for each source

function fetch_burcat(options::Dict=Dict())
    # Implementation details for retrieving data from Burcat database
    url = get(options, "url", "http://garfield.chem.elte.hu/Burcat/THERM.DAT")
    cache_dir = get(options, "cache_dir", joinpath("data", "cache", "burcat"))
    local_file = joinpath(cache_dir, "BURCAT.THR")
    
    # Create directory if needed
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    
    # Download file
    try
        @info "Downloading Burcat database from $(url)"
        response = HTTP.get(url)
        if response.status == 200
            open(local_file, "w") do io
                write(io, response.body)
            end
            @info "Downloaded Burcat database to $(local_file)"
        else
            error("Failed to download Burcat database: HTTP status $(response.status)")
        end
    catch e
        @warn "Error downloading Burcat database: $(e)"
        if !isfile(local_file)
            error("No local Burcat database file available")
        end
    end
    
    # Basic structure - the full implementation would parse the Burcat data
    return Dict(
        "polynomials" => Dict(),  # Will be filled by the actual parser
        "source_path" => local_file,
        "version" => "2023",
        "cas_registry" => Dict()
    )
end

function fetch_nasa(options::Dict=Dict())
    # Implementation details for NASA CEA database
    url = get(options, "url", "https://cearun.grc.nasa.gov/cea/thermo.inp")
    cache_dir = get(options, "cache_dir", joinpath("data", "cache", "nasa-cea"))
    local_file = joinpath(cache_dir, "thermo.inp")
    
    # Create directory if needed
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    
    # Download file
    try
        @info "Downloading NASA CEA database from $(url)"
        response = HTTP.get(url)
        if response.status == 200
            open(local_file, "w") do io
                write(io, response.body)
            end
            @info "Downloaded NASA CEA database to $(local_file)"
        else
            error("Failed to download NASA CEA database: HTTP status $(response.status)")
        end
    catch e
        @warn "Error downloading NASA CEA database: $(e)"
        if !isfile(local_file)
            error("No local NASA CEA database file available")
        end
    end
    
    # Basic structure - the full implementation would parse the NASA data
    return Dict(
        "polynomials" => Dict(),  # Will be filled by the actual parser
        "source_path" => local_file,
        "version" => "2023",
        "cas_registry" => Dict()
    )
end

function fetch_janaf(options::Dict=Dict())
    # Basic structure for JANAF tables
    return Dict(
        "polynomials" => Dict(),
        "version" => "2023"
    )
end

function fetch_nist(options::Dict=Dict())
    # Basic structure for NIST data
    return Dict(
        "polynomials" => Dict(),
        "version" => "2023"
    )
end

function fetch_atct(options::Dict=Dict())
    # Basic structure for ATcT data
    return Dict(
        "polynomials" => Dict(),
        "version" => "2023"
    )
end

function fetch_tde(options::Dict=Dict())
    # Basic structure for TDE data
    return Dict(
        "polynomials" => Dict(),
        "version" => "2023"
    )
end

function fetch_thermoml(options::Dict=Dict())
    # Basic structure for ThermoML data
    return Dict(
        "polynomials" => Dict(),
        "version" => "2023"
    )
end

function fetch_gri_mech(options::Dict=Dict())
    # Basic structure for GRI-Mech 3.0 data
    url = get(options, "url", "http://combustion.berkeley.edu/gri-mech/version30/files30/therm30.dat")
    cache_dir = get(options, "cache_dir", joinpath("data", "cache", "gri-mech"))
    local_file = joinpath(cache_dir, "therm.dat")
    
    # Create directory if needed
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    
    # Download file
    try
        @info "Downloading GRI-Mech database from $(url)"
        response = HTTP.get(url)
        if response.status == 200
            open(local_file, "w") do io
                write(io, response.body)
            end
            @info "Downloaded GRI-Mech database to $(local_file)"
        else
            error("Failed to download GRI-Mech database: HTTP status $(response.status)")
        end
    catch e
        @warn "Error downloading GRI-Mech database: $(e)"
        if !isfile(local_file)
            error("No local GRI-Mech database file available")
        end
    end
    
    return Dict(
        "polynomials" => Dict(),  # Will be filled by the actual parser
        "source_path" => local_file,
        "version" => "3.0"
    )
end

"""
    fetch_all_sources(config::Dict)

Fetch all data sources specified in the configuration.
"""
function fetch_all_sources(config::Dict)
    cache_dir = get(config, "general", Dict())["cache_directory"]
    mkpath(cache_dir)
    
    source_paths = Dict{String, Dict}()
    
    for source in get(config, "data_sources", [])
        if get(source, "enabled", true)
            name = source["name"]
            
            try
                options = Dict(
                    "cache_dir" => joinpath(cache_dir, lowercase(name)),
                    "url" => get(source, "url", ""),
                    "refresh_interval" => get(source, "refresh_interval_days", 7)
                )
                
                source_data = fetch_data_source(name, options)
                source_paths[name] = source_data
                @info "Fetched source: $name"
            catch e
                @warn "Failed to fetch source: $name. Error: $e"
            end
        else
            @info "Skipping disabled source: $(source["name"])"
        end
    end
    
    return source_paths
end