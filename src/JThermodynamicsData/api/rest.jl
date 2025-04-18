"""
REST API for the JThermodynamics package.
This module requires HTTP.jl to be loaded.
"""

"""
    start_rest_api(conn::DuckDB.DB, config::Dict; host::String="0.0.0.0", port::Int=8080)

Start a REST API server for the JThermodynamics package.
"""
function start_rest_api(conn::DuckDB.DB, config::Dict; host::String="0.0.0.0", port::Int=0)
    # Use port from config if not provided
    if port == 0
        port = config["api"]["port"]
    end
    
    # Define routes
    router = HTTP.Router()
    
    # Species endpoint
    HTTP.register!(router, "GET", "/api/species", req -> handle_list_species(req, conn))
    HTTP.register!(router, "GET", "/api/species/:name", req -> handle_get_species(req, conn, config))
    
    # Properties endpoint
    HTTP.register!(router, "GET", "/api/properties/:name", req -> handle_get_properties(req, conn, config))
    
    # Databases endpoint
    HTTP.register!(router, "GET", "/api/databases", req -> handle_list_databases(req, conn))
    
    # CORS middleware
    allowed_origins = config["api"]["allowed_origins"]
    
    cors_middleware = req -> begin
        origin = HTTP.header(req, "Origin", "*")
        
        if "*" in allowed_origins || origin in allowed_origins
            return req
        else
            return HTTP.Response(403, "Origin not allowed")
        end
    end
    
    # Rate limiting middleware
    rate_limit = config["api"]["rate_limit"]
    rate_limiter = RateLimiter(rate_limit)
    
    rate_limit_middleware = req -> begin
        client_ip = HTTP.header(req, "X-Forwarded-For", HTTP.header(req, "X-Real-IP", "127.0.0.1"))
        
        if !check_rate_limit(rate_limiter, client_ip)
            return HTTP.Response(429, "Rate limit exceeded")
        else
            return req
        end
    end
    
    # Logging middleware
    logging_middleware = req -> begin
        @info "$(now()) - $(req.method) $(req.target)"
        return req
    end
    
    # Chain middlewares
    middleware_stack = [logging_middleware, cors_middleware, rate_limit_middleware]
    
    # Start server with middlewares
    @info "Starting REST API server on $host:$port"
    
    HTTP.serve(router, host, port) do req
        # Apply middlewares
        for middleware in middleware_stack
            req = middleware(req)
            if isa(req, HTTP.Response)
                return req
            end
        end
        
        # Add CORS headers to all responses
        try
            resp = HTTP.handle(router, req)
            if "*" in allowed_origins
                HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
            else
                origin = HTTP.header(req, "Origin", "")
                if origin in allowed_origins
                    HTTP.setheader(resp, "Access-Control-Allow-Origin" => origin)
                end
            end
            
            HTTP.setheader(resp, "Access-Control-Allow-Methods" => "GET, OPTIONS")
            HTTP.setheader(resp, "Access-Control-Allow-Headers" => "Content-Type")
            
            return resp
        catch e
            @error "Error handling request: $e"
            return HTTP.Response(500, "Internal server error")
        end
    end
end

"""
    handle_list_species(req::HTTP.Request, conn::DuckDB.DB)

Handle GET /api/species request.
"""
function handle_list_species(req::HTTP.Request, conn::DuckDB.DB)
    # Parse query parameters
    params = HTTP.queryparams(HTTP.URI(req.target))
    
    limit = 100
    if haskey(params, "limit")
        try
            limit = parse(Int, params["limit"])
            limit = min(max(limit, 1), 1000)  # Limit between 1 and 1000
        catch
            # Ignore parse errors
        end
    end
    
    search = get(params, "search", "")
    
    # Query species
    species = list_available_species(conn, limit=limit, search_term=search)
    
    # Convert to JSON
    result = Dict(
        "species" => [Dict(pairs(row)) for row in eachrow(species)],
        "count" => size(species, 1),
        "limit" => limit
    )
    
    return HTTP.Response(200, JSON.json(result))
end

"""
    handle_get_species(req::HTTP.Request, conn::DuckDB.DB, config::Dict)

Handle GET /api/species/:name request.
"""
function handle_get_species(req::HTTP.Request, conn::DuckDB.DB, config::Dict)
    # Extract species name from URL
    m = match(r"/api/species/([^/]+)", req.target)
    
    if m === nothing
        return HTTP.Response(400, "Invalid species name")
    end
    
    species_name = HTTP.unescapeuri(m.captures[1])
    
    # Query species
    try
        species = query_species(conn, species_name, config)
        return HTTP.Response(200, JSON.json(species))
    catch e
        if isa(e, ErrorException) && startswith(e.msg, "Species not found")
            return HTTP.Response(404, "Species not found: $species_name")
        else
            @error "Error querying species: $e"
            return HTTP.Response(500, "Internal server error")
        end
    end
end

"""
    handle_get_properties(req::HTTP.Request, conn::DuckDB.DB, config::Dict)

Handle GET /api/properties/:name request.
"""
function handle_get_properties(req::HTTP.Request, conn::DuckDB.DB, config::Dict)
    # Extract species name from URL
    m = match(r"/api/properties/([^/]+)", req.target)
    
    if m === nothing
        return HTTP.Response(400, "Invalid species name")
    end
    
    species_name = HTTP.unescapeuri(m.captures[1])
    
    # Parse query parameters
    params = HTTP.queryparams(HTTP.URI(req.target))
    
    # Get temperature
    temperature = 298.15
    if haskey(params, "temperature")
        try
            temperature = parse(Float64, params["temperature"])
        catch
            return HTTP.Response(400, "Invalid temperature parameter")
        end
    end
    
    # Check temperature range
    temp_range = config["general"]["temperature_range"]
    if temperature < temp_range[1] || temperature > temp_range[2]
        return HTTP.Response(400, "Temperature out of range ($(temp_range[1])-$(temp_range[2]) K)")
    end
    
    # Calculate range of temperatures
    calculate_range = false
    if haskey(params, "range") && params["range"] == "true"
        calculate_range = true
    end
    
    # Parse range parameters
    temp_min = temperature
    temp_max = temperature
    temp_step = 10.0
    
    if calculate_range
        if haskey(params, "temp_min")
            try
                temp_min = parse(Float64, params["temp_min"])
            catch
                return HTTP.Response(400, "Invalid temp_min parameter")
            end
        end
        
        if haskey(params, "temp_max")
            try
                temp_max = parse(Float64, params["temp_max"])
            catch
                return HTTP.Response(400, "Invalid temp_max parameter")
            end
        end
        
        if haskey(params, "temp_step")
            try
                temp_step = parse(Float64, params["temp_step"])
            catch
                return HTTP.Response(400, "Invalid temp_step parameter")
            end
        end
        
        # Validate range
        if temp_min > temp_max
            return HTTP.Response(400, "temp_min must be less than or equal to temp_max")
        end
        
        if temp_min < temp_range[1] || temp_max > temp_range[2]
            return HTTP.Response(400, "Temperature range out of bounds ($(temp_range[1])-$(temp_range[2]) K)")
        end
        
        if temp_step <= 0 || temp_step > (temp_max - temp_min)
            return HTTP.Response(400, "Invalid temp_step parameter")
        end
    end
    
    # Query properties
    try
        if calculate_range
            result = calculate_properties(conn, species_name, [temp_min, temp_max], config, step=temp_step)
        else
            result = query_properties(conn, species_name, temperature, config)
        end
        
        return HTTP.Response(200, JSON.json(result))
    catch e
        if isa(e, ErrorException) && startswith(e.msg, "Species not found")
            return HTTP.Response(404, "Species not found: $species_name")
        elseif isa(e, ErrorException) && startswith(e.msg, "No data found")
            return HTTP.Response(404, e.msg)
        else
            @error "Error querying properties: $e"
            return HTTP.Response(500, "Internal server error")
        end
    end
end

"""
    handle_list_databases(req::HTTP.Request, conn::DuckDB.DB)

Handle GET /api/databases request.
"""
function handle_list_databases(req::HTTP.Request, conn::DuckDB.DB)
    # Query databases
    databases = list_available_databases(conn)
    
    # Convert to JSON
    result = Dict(
        "databases" => [Dict(pairs(row)) for row in eachrow(databases)],
        "count" => size(databases, 1)
    )
    
    return HTTP.Response(200, JSON.json(result))
end

"""
Rate limiter for API requests.
"""
mutable struct RateLimiter
    limit::Int
    window::Int
    clients::Dict{String, Vector{Float64}}
    
    function RateLimiter(limit::Int, window::Int=60)
        new(limit, window, Dict{String, Vector{Float64}}())
    end
end

"""
    check_rate_limit(limiter::RateLimiter, client_ip::String)

Check if a client has exceeded their rate limit.
"""
function check_rate_limit(limiter::RateLimiter, client_ip::String)
    now_time = time()
    
    # Initialize client if not exists
    if !haskey(limiter.clients, client_ip)
        limiter.clients[client_ip] = Float64[]
    end
    
    # Remove old timestamps
    window_start = now_time - limiter.window
    limiter.clients[client_ip] = filter(t -> t >= window_start, limiter.clients[client_ip])
    
    # Check if under limit
    if length(limiter.clients[client_ip]) < limiter.limit
        push!(limiter.clients[client_ip], now_time)
        return true
    else
        return false
    end
end