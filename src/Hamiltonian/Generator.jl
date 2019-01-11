"""
Generator submodule
"""
module Generator

using LinearAlgebra

using ..Lattice
using ..Spec


"""
    hopping_inplace
"""
function hopping_inplace(uc ::UnitCell{O}, hop::HoppingDiagonal{R}) where {O, R<:Real}
    ndim = dimension(uc)
    norb = numorbital(uc)
    v = hop.amplitude
    i = hop.i

    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        @assert(size(momentum) == (ndim,))
        @assert(size(out) == (norb, norb))
        out[i,i] += v
        return out
    end
end


"""
    hopping_inplace
"""
function hopping_inplace(uc ::UnitCell{O}, hop::HoppingOffdiagonal{C}) where {O, C<:Number}
    ndim = dimension(uc)
    norb = numorbital(uc)
    v = hop.amplitude
    (i, j) = (hop.i, hop.j)
    ri = fract2carte(uc, getorbitalcoord(uc, hop.i) + hop.Ri)
    rj = fract2carte(uc, getorbitalcoord(uc, hop.j) + hop.Rj)
    rij = rj - ri

    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        @assert(size(momentum) == (ndim,))
        @assert(size(out) == (norb, norb))
        phase = cis(dot(momentum, rij))
        out[i,j] += v * phase
        out[j,i] += conj(v * phase)
        return out
    end
end


#=
function hopping_inplace(uc ::UnitCell{O},
                         hops ::AbstractVector{Hopping}) where {O}
    ndim = dimension(uc)
    norb = numorbital(uc)
    funcs = [hopping_inplace(uc, hop) for hop in hops]

    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        for func in funcs
            func(momentum, out)
        end
        return out
    end
end
=#

"""
    hopping_inplace
"""
function hopping_inplace(uc ::UnitCell{O},
                         hops ::AbstractVector{Hopping}) where {O}
    hops_diag = HoppingDiagonal[]
    hops_offdiag = HoppingOffdiagonal[]
    for hop in hops
        if isa(hop, HoppingDiagonal)
            push!(hops_diag, hop)
        elseif isa(hop, HoppingOffdiagonal)
            push!(hops_offdiag, hop)
        else
            assert(false, "Hopping should be either diagonal or offdiagonal")
        end
    end
    func_diag = hopping_inplace(uc, hops_diag)
    func_offdiag = hopping_inplace(uc, hops_offdiag)
    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        func_diag(momentum, out)
        func_offdiag(momentum, out)
    end
end

"""
    hopping_inplace
"""
function hopping_inplace(uc ::UnitCell{O},
                         hops ::AbstractVector{HoppingDiagonal}) where {O}
    ndim = dimension(uc)
    norb = numorbital(uc)
    hops_diag = Tuple{Float64, Int}[]

    for hop in hops
        v = hop.amplitude
        i = hop.i
        push!(hops_diag, (v, i))
    end

    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        @assert(size(momentum) == (ndim,))
        @assert(size(out) == (norb, norb))
        for (v, i) in hops_diag
            out[i,i] += v
        end
        return out
    end
end


"""
    hopping_inplace
"""
function hopping_inplace(uc ::UnitCell{O},
                         hops ::AbstractVector{HoppingOffdiagonal}) where {O}
    ndim = dimension(uc)
    norb = numorbital(uc)
    hops_offdiag = Tuple{ComplexF64, Int, Int, Vector{Float64}}[]

    for hop in hops
        v = hop.amplitude
        (i, j) = (hop.i, hop.j)
        ri = fract2carte(uc, getorbitalcoord(uc, hop.i) + hop.Ri)
        rj = fract2carte(uc, getorbitalcoord(uc, hop.j) + hop.Rj)
        rij = rj - ri
        push!(hops_offdiag, (v, i, j, rij))
    end

    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        @assert(size(momentum) == (ndim,))
        @assert(size(out) == (norb, norb))
        for (v, i, j, rij) in hops_offdiag
            phase = cis(dot(momentum, rij))
            v_phase = v * phase
            out[i,j] += v_phase
            out[j,i] += conj(v_phase)
        end
        return out
    end
end


"""
    hopping_inplace
"""
function hopping_inplace(uc ::UnitCell{O},
                         the_hops_diag ::AbstractVector{HoppingDiagonal},
                         the_hops_offdiag ::AbstractVector{HoppingOffdiagonal}) where {O}
    ndim = dimension(uc)
    norb = numorbital(uc)
    hops_diag = Tuple{Float64, Int}[]
    hops_offdiag = Tuple{ComplexF64, Int, Int, Vector{Float64}}[]

    for hop in the_hops_diag
        v = hop.amplitude
        i = hop.i
        push!(hops_diag, (v, i))
    end

    for hop in the_hops_offdiag
        v = hop.amplitude
        (i, j) = (hop.i, hop.j)
        ri = fract2carte(uc, getorbitalcoord(uc, hop.i) + hop.Ri)
        rj = fract2carte(uc, getorbitalcoord(uc, hop.j) + hop.Rj)
        rij = rj - ri
        push!(hops_offdiag, (v, i, j, rij))
    end

    function(momentum ::AbstractVector{Float64}, out::AbstractArray{ComplexF64, 2})
        @assert(size(momentum) == (ndim,))
        @assert(size(out) == (norb, norb))
        for (v, i) in hops_diag
            out[i,i] += v
        end
        for (v, i, j, rij) in hops_offdiag
            phase = cis(dot(momentum, rij))
            v_phase = v * phase
            out[i,j] += v_phase
            out[j,i] += conj(v_phase)
        end
        return out
    end
end


"""
    hopping_inplace
"""
function hopping_inplace(hamiltonian ::FullHamiltonian{O}) where {O}
    return hopping_inplace(hamiltonian.unitcell,
                           hamiltonian.hoppings_diagonal,
                           hamiltonian.hoppings_offdiagonal)
end

end
