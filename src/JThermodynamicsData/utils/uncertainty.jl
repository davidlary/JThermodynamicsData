"""
Functions for handling uncertainty propagation in thermodynamic calculations.
"""

using Measurements
using MonteCarloMeasurements

"""
    apply_uncertainty(value::Float64, uncertainty::Float64, method::UncertaintyMethod)

Apply uncertainty to a value using the specified method.
"""
function apply_uncertainty(value::Float64, uncertainty::Float64, method::UncertaintyMethod)
    if method == MEASUREMENTS
        return value ± uncertainty
    elseif method == MONTE_CARLO
        return Particles(1000, value, uncertainty)
    elseif method == INTERVAL
        # Interval arithmetic (simplified)
        return value  # Return just the central value for now
    else
        return value  # No uncertainty
    end
end

"""
    extract_uncertainty(value::Union{Float64, Measurement{Float64}, Particles{Float64}})

Extract the uncertainty from a value.
"""
function extract_uncertainty(value::Union{Float64, Measurement{Float64}, Particles{Float64}})
    if isa(value, Measurement)
        return Measurements.uncertainty(value)
    elseif isa(value, Particles)
        return std(value)
    else
        return 0.0  # No uncertainty
    end
end

"""
    central_value(value::Union{Float64, Measurement{Float64}, Particles{Float64}})

Extract the central value from a value with uncertainty.
"""
function central_value(value::Union{Float64, Measurement{Float64}, Particles{Float64}})
    if isa(value, Measurement)
        return Measurements.value(value)
    elseif isa(value, Particles)
        return mean(value)
    else
        return value
    end
end

"""
    combine_uncertainties(values::Vector{Float64}, uncertainties::Vector{Float64}, 
                        weights::Vector{Float64}=Float64[])

Combine multiple values with uncertainties using weighted averaging.
"""
function combine_uncertainties(values::Vector{Float64}, uncertainties::Vector{Float64}, 
                             weights::Vector{Float64}=Float64[])
    if length(values) != length(uncertainties)
        error("Values and uncertainties must have the same length")
    end
    
    if isempty(values)
        return 0.0, 0.0
    end
    
    # If weights are not provided, use inverse variance weighting
    if isempty(weights)
        # Calculate weights as 1/variance (with safeguards for zero uncertainty)
        weights = [u > 0 ? 1/(u^2) : 1.0 for u in uncertainties]
    end
    
    # Normalize weights
    total_weight = sum(weights)
    if total_weight == 0
        # Equal weighting if all weights are zero
        weights = ones(length(values)) / length(values)
    else
        weights = weights / total_weight
    end
    
    # Calculate weighted average
    avg = sum(values .* weights)
    
    # Calculate combined uncertainty
    if length(values) == 1
        # Single value, just return its uncertainty
        unc = uncertainties[1]
    else
        # Use weighted standard deviation for combined uncertainty
        # First calculate weighted variance
        variance = sum(weights .* (values .- avg).^2)
        
        # Scale factor for sample variance
        scale = length(values) / (length(values) - 1)
        
        # Combined uncertainty is square root of weighted variance
        unc = sqrt(variance * scale)
        
        # Add in uncertainty from the individual measurements
        # sqrt(sum of squared weighted uncertainties)
        measurement_unc = sqrt(sum((weights .* uncertainties).^2))
        
        # Total uncertainty is root sum of squares of the two components
        unc = sqrt(unc^2 + measurement_unc^2)
    end
    
    return avg, unc
end

"""
    propagate_uncertainty_cp(coeffs::Vector{Float64}, coeffs_uncertainty::Vector{Float64}, 
                          temperature::Float64, method::UncertaintyMethod)

Propagate uncertainty in heat capacity calculation.
"""
function propagate_uncertainty_cp(coeffs::Vector{Float64}, coeffs_uncertainty::Vector{Float64}, 
                               temperature::Float64, method::UncertaintyMethod)
    if length(coeffs) != length(coeffs_uncertainty)
        error("Coefficients and uncertainties must have the same length")
    end
    
    if method == MEASUREMENTS
        # Use Measurements.jl for linear uncertainty propagation
        coeffs_with_uncertainty = [coeffs[i] ± coeffs_uncertainty[i] for i in 1:length(coeffs)]
        
        if length(coeffs) == 7  # NASA 7
            return coeffs_with_uncertainty[1] + 
                   coeffs_with_uncertainty[2] * temperature + 
                   coeffs_with_uncertainty[3] * temperature^2 + 
                   coeffs_with_uncertainty[4] * temperature^3 + 
                   coeffs_with_uncertainty[5] * temperature^4
        elseif length(coeffs) == 9  # NASA 9
            return coeffs_with_uncertainty[1] * temperature^(-2) + 
                   coeffs_with_uncertainty[2] * temperature^(-1) + 
                   coeffs_with_uncertainty[3] + 
                   coeffs_with_uncertainty[4] * temperature + 
                   coeffs_with_uncertainty[5] * temperature^2 + 
                   coeffs_with_uncertainty[6] * temperature^3 + 
                   coeffs_with_uncertainty[7] * temperature^4
        else
            error("Unsupported coefficient length")
        end
    elseif method == MONTE_CARLO
        # Use MonteCarloMeasurements.jl for Monte Carlo uncertainty propagation
        coeffs_particles = [Particles(1000, coeffs[i], coeffs_uncertainty[i]) for i in 1:length(coeffs)]
        
        if length(coeffs) == 7  # NASA 7
            return coeffs_particles[1] + 
                   coeffs_particles[2] * temperature + 
                   coeffs_particles[3] * temperature^2 + 
                   coeffs_particles[4] * temperature^3 + 
                   coeffs_particles[5] * temperature^4
        elseif length(coeffs) == 9  # NASA 9
            return coeffs_particles[1] * temperature^(-2) + 
                   coeffs_particles[2] * temperature^(-1) + 
                   coeffs_particles[3] + 
                   coeffs_particles[4] * temperature + 
                   coeffs_particles[5] * temperature^2 + 
                   coeffs_particles[6] * temperature^3 + 
                   coeffs_particles[7] * temperature^4
        else
            error("Unsupported coefficient length")
        end
    else  # INTERVAL or NONE
        # Calculate central value
        if length(coeffs) == 7  # NASA 7
            cp = coeffs[1] + 
                 coeffs[2] * temperature + 
                 coeffs[3] * temperature^2 + 
                 coeffs[4] * temperature^3 + 
                 coeffs[5] * temperature^4
        elseif length(coeffs) == 9  # NASA 9
            cp = coeffs[1] * temperature^(-2) + 
                 coeffs[2] * temperature^(-1) + 
                 coeffs[3] + 
                 coeffs[4] * temperature + 
                 coeffs[5] * temperature^2 + 
                 coeffs[6] * temperature^3 + 
                 coeffs[7] * temperature^4
        else
            error("Unsupported coefficient length")
        end
        
        # Calculate uncertainty using sensitivity coefficients
        unc_sq = 0.0
        
        if length(coeffs) == 7  # NASA 7
            sensitivities = [
                1.0,  # a1
                temperature,  # a2
                temperature^2,  # a3
                temperature^3,  # a4
                temperature^4   # a5
            ]
            
            for i in 1:5
                unc_sq += (sensitivities[i] * coeffs_uncertainty[i])^2
            end
        elseif length(coeffs) == 9  # NASA 9
            sensitivities = [
                temperature^(-2),  # a1
                temperature^(-1),  # a2
                1.0,              # a3
                temperature,      # a4
                temperature^2,    # a5
                temperature^3,    # a6
                temperature^4     # a7
            ]
            
            for i in 1:7
                unc_sq += (sensitivities[i] * coeffs_uncertainty[i])^2
            end
        end
        
        # Return central value and uncertainty
        return cp, sqrt(unc_sq)
    end
end

"""
    propagate_uncertainty_h(coeffs::Vector{Float64}, coeffs_uncertainty::Vector{Float64}, 
                         temperature::Float64, method::UncertaintyMethod)

Propagate uncertainty in enthalpy calculation.
"""
function propagate_uncertainty_h(coeffs::Vector{Float64}, coeffs_uncertainty::Vector{Float64}, 
                              temperature::Float64, method::UncertaintyMethod)
    if length(coeffs) != length(coeffs_uncertainty)
        error("Coefficients and uncertainties must have the same length")
    end
    
    if method == MEASUREMENTS
        # Use Measurements.jl for linear uncertainty propagation
        coeffs_with_uncertainty = [coeffs[i] ± coeffs_uncertainty[i] for i in 1:length(coeffs)]
        
        if length(coeffs) == 7  # NASA 7
            return coeffs_with_uncertainty[1] + 
                   coeffs_with_uncertainty[2] * temperature / 2 + 
                   coeffs_with_uncertainty[3] * temperature^2 / 3 + 
                   coeffs_with_uncertainty[4] * temperature^3 / 4 + 
                   coeffs_with_uncertainty[5] * temperature^4 / 5 + 
                   coeffs_with_uncertainty[6] / temperature
        elseif length(coeffs) == 9  # NASA 9
            return -coeffs_with_uncertainty[1] * temperature^(-2) + 
                   coeffs_with_uncertainty[2] * log(temperature) / temperature + 
                   coeffs_with_uncertainty[3] + 
                   coeffs_with_uncertainty[4] * temperature / 2 + 
                   coeffs_with_uncertainty[5] * temperature^2 / 3 + 
                   coeffs_with_uncertainty[6] * temperature^3 / 4 + 
                   coeffs_with_uncertainty[7] * temperature^4 / 5 + 
                   coeffs_with_uncertainty[8] / temperature
        else
            error("Unsupported coefficient length")
        end
    elseif method == MONTE_CARLO
        # Use MonteCarloMeasurements.jl for Monte Carlo uncertainty propagation
        coeffs_particles = [Particles(1000, coeffs[i], coeffs_uncertainty[i]) for i in 1:length(coeffs)]
        
        if length(coeffs) == 7  # NASA 7
            return coeffs_particles[1] + 
                   coeffs_particles[2] * temperature / 2 + 
                   coeffs_particles[3] * temperature^2 / 3 + 
                   coeffs_particles[4] * temperature^3 / 4 + 
                   coeffs_particles[5] * temperature^4 / 5 + 
                   coeffs_particles[6] / temperature
        elseif length(coeffs) == 9  # NASA 9
            return -coeffs_particles[1] * temperature^(-2) + 
                   coeffs_particles[2] * log(temperature) / temperature + 
                   coeffs_particles[3] + 
                   coeffs_particles[4] * temperature / 2 + 
                   coeffs_particles[5] * temperature^2 / 3 + 
                   coeffs_particles[6] * temperature^3 / 4 + 
                   coeffs_particles[7] * temperature^4 / 5 + 
                   coeffs_particles[8] / temperature
        else
            error("Unsupported coefficient length")
        end
    else  # INTERVAL or NONE
        # Calculate central value
        if length(coeffs) == 7  # NASA 7
            h = coeffs[1] + 
                coeffs[2] * temperature / 2 + 
                coeffs[3] * temperature^2 / 3 + 
                coeffs[4] * temperature^3 / 4 + 
                coeffs[5] * temperature^4 / 5 + 
                coeffs[6] / temperature
        elseif length(coeffs) == 9  # NASA 9
            h = -coeffs[1] * temperature^(-2) + 
                coeffs[2] * log(temperature) / temperature + 
                coeffs[3] + 
                coeffs[4] * temperature / 2 + 
                coeffs[5] * temperature^2 / 3 + 
                coeffs[6] * temperature^3 / 4 + 
                coeffs[7] * temperature^4 / 5 + 
                coeffs[8] / temperature
        else
            error("Unsupported coefficient length")
        end
        
        # Calculate uncertainty using sensitivity coefficients
        unc_sq = 0.0
        
        if length(coeffs) == 7  # NASA 7
            sensitivities = [
                1.0,                  # a1
                temperature / 2,      # a2
                temperature^2 / 3,    # a3
                temperature^3 / 4,    # a4
                temperature^4 / 5,    # a5
                1 / temperature       # a6
            ]
            
            for i in 1:6
                unc_sq += (sensitivities[i] * coeffs_uncertainty[i])^2
            end
        elseif length(coeffs) == 9  # NASA 9
            sensitivities = [
                -temperature^(-2),                # a1
                log(temperature) / temperature,   # a2
                1.0,                             # a3
                temperature / 2,                 # a4
                temperature^2 / 3,               # a5
                temperature^3 / 4,               # a6
                temperature^4 / 5,               # a7
                1 / temperature                  # a8
            ]
            
            for i in 1:8
                unc_sq += (sensitivities[i] * coeffs_uncertainty[i])^2
            end
        end
        
        # Return central value and uncertainty
        return h, sqrt(unc_sq)
    end
end

"""
    propagate_uncertainty_s(coeffs::Vector{Float64}, coeffs_uncertainty::Vector{Float64}, 
                         temperature::Float64, method::UncertaintyMethod)

Propagate uncertainty in entropy calculation.
"""
function propagate_uncertainty_s(coeffs::Vector{Float64}, coeffs_uncertainty::Vector{Float64}, 
                              temperature::Float64, method::UncertaintyMethod)
    if length(coeffs) != length(coeffs_uncertainty)
        error("Coefficients and uncertainties must have the same length")
    end
    
    if method == MEASUREMENTS
        # Use Measurements.jl for linear uncertainty propagation
        coeffs_with_uncertainty = [coeffs[i] ± coeffs_uncertainty[i] for i in 1:length(coeffs)]
        
        if length(coeffs) == 7  # NASA 7
            return coeffs_with_uncertainty[1] * log(temperature) + 
                   coeffs_with_uncertainty[2] * temperature + 
                   coeffs_with_uncertainty[3] * temperature^2 / 2 + 
                   coeffs_with_uncertainty[4] * temperature^3 / 3 + 
                   coeffs_with_uncertainty[5] * temperature^4 / 4 + 
                   coeffs_with_uncertainty[7]
        elseif length(coeffs) == 9  # NASA 9
            return -coeffs_with_uncertainty[1] * temperature^(-2) / 2 - 
                   coeffs_with_uncertainty[2] * temperature^(-1) + 
                   coeffs_with_uncertainty[3] * log(temperature) + 
                   coeffs_with_uncertainty[4] * temperature + 
                   coeffs_with_uncertainty[5] * temperature^2 / 2 + 
                   coeffs_with_uncertainty[6] * temperature^3 / 3 + 
                   coeffs_with_uncertainty[7] * temperature^4 / 4 + 
                   coeffs_with_uncertainty[9]
        else
            error("Unsupported coefficient length")
        end
    elseif method == MONTE_CARLO
        # Use MonteCarloMeasurements.jl for Monte Carlo uncertainty propagation
        coeffs_particles = [Particles(1000, coeffs[i], coeffs_uncertainty[i]) for i in 1:length(coeffs)]
        
        if length(coeffs) == 7  # NASA 7
            return coeffs_particles[1] * log(temperature) + 
                   coeffs_particles[2] * temperature + 
                   coeffs_particles[3] * temperature^2 / 2 + 
                   coeffs_particles[4] * temperature^3 / 3 + 
                   coeffs_particles[5] * temperature^4 / 4 + 
                   coeffs_particles[7]
        elseif length(coeffs) == 9  # NASA 9
            return -coeffs_particles[1] * temperature^(-2) / 2 - 
                   coeffs_particles[2] * temperature^(-1) + 
                   coeffs_particles[3] * log(temperature) + 
                   coeffs_particles[4] * temperature + 
                   coeffs_particles[5] * temperature^2 / 2 + 
                   coeffs_particles[6] * temperature^3 / 3 + 
                   coeffs_particles[7] * temperature^4 / 4 + 
                   coeffs_particles[9]
        else
            error("Unsupported coefficient length")
        end
    else  # INTERVAL or NONE
        # Calculate central value
        if length(coeffs) == 7  # NASA 7
            s = coeffs[1] * log(temperature) + 
                coeffs[2] * temperature + 
                coeffs[3] * temperature^2 / 2 + 
                coeffs[4] * temperature^3 / 3 + 
                coeffs[5] * temperature^4 / 4 + 
                coeffs[7]
        elseif length(coeffs) == 9  # NASA 9
            s = -coeffs[1] * temperature^(-2) / 2 - 
                coeffs[2] * temperature^(-1) + 
                coeffs[3] * log(temperature) + 
                coeffs[4] * temperature + 
                coeffs[5] * temperature^2 / 2 + 
                coeffs[6] * temperature^3 / 3 + 
                coeffs[7] * temperature^4 / 4 + 
                coeffs[9]
        else
            error("Unsupported coefficient length")
        end
        
        # Calculate uncertainty using sensitivity coefficients
        unc_sq = 0.0
        
        if length(coeffs) == 7  # NASA 7
            sensitivities = [
                log(temperature),      # a1
                temperature,          # a2
                temperature^2 / 2,    # a3
                temperature^3 / 3,    # a4
                temperature^4 / 4,    # a5
                0.0,                 # a6 (not used in entropy)
                1.0                  # a7
            ]
            
            for i in [1, 2, 3, 4, 5, 7]  # Skip a6
                unc_sq += (sensitivities[i] * coeffs_uncertainty[i])^2
            end
        elseif length(coeffs) == 9  # NASA 9
            sensitivities = [
                -temperature^(-2) / 2,   # a1
                -temperature^(-1),      # a2
                log(temperature),       # a3
                temperature,           # a4
                temperature^2 / 2,     # a5
                temperature^3 / 3,     # a6
                temperature^4 / 4,     # a7
                0.0,                  # a8 (not used in entropy)
                1.0                   # a9
            ]
            
            for i in [1, 2, 3, 4, 5, 6, 7, 9]  # Skip a8
                unc_sq += (sensitivities[i] * coeffs_uncertainty[i])^2
            end
        end
        
        # Return central value and uncertainty
        return s, sqrt(unc_sq)
    end
end