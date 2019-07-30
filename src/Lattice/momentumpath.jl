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


#function momentumpath(mp ::MomentumPath, segment_length ::Integer)
#    ks = linpath(mp.anchors, segment_length)
#    dks = ks[2:end, :] - ks[1:end-1, :]
#    dls = zeros(Float64, size(dks)[1])
#    for i in 1:size(dks)[1]
#        dls[i] = norm(dks[i, :])
#    end
#    ls = Float64[0.0, cumsum(dls)...]
#    ticks = Float64[]
#    ticklabels = String[]
#    for i in 1:length(mp.anchors)
#        push!(ticks, ls[segment_length*(i-1)+1])
#        push!(ticklabels, mp.anchornames[i])
#    end
#    
#    Dict(
#        :momentum=> ks,
#        :momentum_distance => ls,
#        :ticks => ticks,
#        :ticklabels => ticklabels,
#    )
#end
