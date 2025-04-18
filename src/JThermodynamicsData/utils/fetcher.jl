"""
Functions for fetching thermodynamic data from online sources.
"""

"""
    fetch_data_source(source::Dict, cache_dir::String)

Fetch data from a source specified in the configuration.
"""
function fetch_data_source(source::Dict, cache_dir::String)
    name = source["name"]
    format = source["format"]
    url = get(source, "url", "")
    
    # Skip if no URL provided
    if isempty(url)
        @info "No URL provided for source: $name"
        return get(source, "path", "")
    end
    
    # Create cache directory
    source_cache_dir = joinpath(cache_dir, lowercase(name))
    mkpath(source_cache_dir)
    
    # Generate output path based on format
    filename = ""
    
    if format == "chemkin"
        filename = "therm.dat"
    elseif format == "nasa7-9"
        filename = "thermo.inp"
    elseif format == "burcat"
        filename = "BURCAT.THR"
    elseif format == "thermoml"
        filename = "thermoml.xml"
    elseif format == "atct"
        filename = "ATcT.nist"
    else
        # Generic filename
        filename = "data.txt"
    end
    
    output_path = joinpath(source_cache_dir, filename)
    
    # Check if we need to download
    need_download = true
    
    if isfile(output_path)
        # Check file age
        file_age = time() - mtime(output_path)
        refresh_interval = get(source, "refresh_interval_days", 90) * 24 * 60 * 60  # Convert days to seconds
        
        if file_age < refresh_interval
            @info "Cache is fresh for source: $name (age: $(round(file_age / (24 * 60 * 60), digits=1)) days)"
            need_download = false
        end
    end
    
    if need_download
        try
            @info "Downloading data for source: $name"
            download_data_source(url, output_path, format)
            @info "Downloaded data to: $output_path"
        catch e
            @warn "Failed to download data for source: $name. Error: $e"
            
            if isfile(output_path)
                @info "Using existing cache file: $output_path"
                return output_path
            else
                error("No cache file available for source: $name")
            end
        end
    end
    
    return output_path
end

"""
    download_data_source(url::String, output_path::String, format::String)

Download data from a URL and save it to the specified path.
"""
function download_data_source(url::String, output_path::String, format::String)
    # Handle different download methods based on the source format
    if format == "burcat"
        # Use specific burcat downloader
        return download_burcat_database(dirname(output_path))
    elseif format == "atct"
        # Use specific ATcT downloader
        return download_atct_data(dirname(output_path))
    else
        # Generic download with HTTP.jl
        try
            response = HTTP.get(url)
            
            if response.status == 200
                open(output_path, "w") do io
                    write(io, response.body)
                end
                return output_path
            else
                error("Failed to download data: HTTP status $(response.status)")
            end
        catch e
            error("Failed to download data: $e")
        end
    end
end

"""
    fetch_all_sources(config::Dict)

Fetch all data sources specified in the configuration.
"""
function fetch_all_sources(config::Dict)
    cache_dir = config["general"]["cache_directory"]
    mkpath(cache_dir)
    
    source_paths = Dict{String, String}()
    
    for source in config["data_sources"]
        if source["enabled"]
            name = source["name"]
            
            try
                source_path = fetch_data_source(source, cache_dir)
                source_paths[name] = source_path
                @info "Fetched source: $name -> $source_path"
            catch e
                @warn "Failed to fetch source: $name. Error: $e"
            end
        else
            @info "Skipping disabled source: $(source["name"])"
        end
    end
    
    return source_paths
end