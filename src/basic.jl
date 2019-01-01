export fermidirac

"""
    fermidirac

Return the Fermi-Dirac distribution function for the given temperature.
For energy whose absolute value is less than `etol`, return 0.5.
"""
function fermidirac(temperature ::T;
                    ttol::T=eps(T),
                    etol::T=sqrt(eps(T))) where {T <:AbstractFloat}
    if ttol < 0
        throw(ArgumentError("ttol should be non-negative"))
    elseif etol < 0
        throw(ArgumentError("etol should be non-negative"))
    elseif temperature < 0
        throw(ArgumentError("temperature should be non-negative"))
    end
    if temperature < ttol
        function(e ::T)
            if e <= -etol
                return T(1)
            elseif e < etol
                return T(0.5)
            else
                return T(0)
            end
        end
    else
        beta = 1.0 / temperature
        function(e ::T)
            return T(1 / (exp(beta * e) + 1))
        end
    end
end


"""
fermidirac

Return the Fermi-Dirac distribution function for the given temperature.
For energy whose absolute value is less than `etol`, return 0.5.
"""
function fermidirac(temperature ::T;
                    ttol::Float64=eps(Float64),
                    etol::Float64=sqrt(eps(Float64))) where {T <:Integer}
    if ttol < 0
        throw(ArgumentError("ttol should be non-negative"))
    elseif etol < 0
        throw(ArgumentError("etol should be non-negative"))
    elseif temperature < 0
        throw(ArgumentError("temperature should be non-negative"))
    end
    if temperature == 0
        function(e ::Float64)
            if e <= -etol
                return 1.0
            elseif e < etol
                return 0.5
            else
                return 0.0
            end
        end
    else
        beta = 1 / Float64(temperature)
        function(e ::Float64)
            return 1 / (exp(beta * e) + 1)
        end
    end
end
