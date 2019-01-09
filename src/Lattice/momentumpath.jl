export linpath
export momentumpath

"""
    momentumpath

Generate a list of momenta

# Arguments
* `anchorpoints`
* (Optional) `nseg` - number of points in each segment
"""
function linpath(anchorpoints ::Vector{Vector{Float64}}; nseg=64)
    n = length(anchorpoints)
    if n <= 1
        throw(ArgumentError("Number of anchor points should be more than one"))
    end

    d = length(anchorpoints[1])
    for i in 2:n
        if length(anchorpoints[i]) != d
            throw(ArgumentError("All anchor points should have the same dimension"))
        end
    end

    out = Vector{Float64}[]
    for i in 1:n-1
        from = anchorpoints[i]
        to = anchorpoints[i+1]
        dvec = (to - from) ./ nseg
        append!(out, from + dvec * 0:(nseg-1))
    end
    push!(out, anchorpoints[end])
    return out
end


"""
    momentumpath

The anchorpoints are given in units of the reciprocal lattice vectors.
"""
function momentumpath(unitcell ::UnitCell,
                      anchorpoints ::Vector{Vector{Float64}})
    real_anchorpoints = [unitcell.reciprocallatticevectors * ap for ap in anchorpoints]
    return linpath(real_anchorpoints)
end
