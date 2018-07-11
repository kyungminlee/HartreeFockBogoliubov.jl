
export isvalidtimereversalmatrix

export istimereversal, iscanonicaltimereversalinvariant, timereversalindexgrid, degenerategroups
export kramerpairup!, kramerpairupblock!

if VERSION < v"0.7-"
    using MicroLogging
end

"""
    isvalidtimereversalmatrix

Test whether the given matrix is a valid unitary matrix for the time reversal operation.

    ```math
    T = U ⋅ K
    ```

    ``U`` must satisfy the two conditions:
    1. ``U U^{\\dagger} = 1`` (from unitarity of ``U``)
    2. ``U = - U^{\\mathsf{T}}`` (from `T^2 = -1`)
"""
function isvalidtimereversalmatrix(
            mat::AbstractMatrix{<:Number};
            tol::Real=sqrt(eps(Float64)))
    (nr, nc) = size(mat)
    if nr != nc
        return false
    elseif ! isapprox(mat + mat.', zeros(mat); atol=tol)
        return false
    elseif ! isapprox(mat * mat', eye(mat); atol=tol)
        return false
    else
        return true
    end
end


function istimereversal(a ::AbstractMatrix{C2}, theta ::AbstractMatrix{C1}) where {C1 <: Number, C2 <: Number}
  @show maximum(abs.( a - conj(theta' * a * theta)))
  return all(isapprox.( a - conj(theta' * a * theta), 0
                          ; atol=sqrt(eps(Float64))))
end



function iscanonicaltimereversalinvariant(
    unitcell::UnitCell,
    timereversalmatrix::AbstractMatrix{<:Number},
    hoppings::AbstractVector{Hopping};
    tol::Real = sqrt(eps(Float64)),
    momentumgridsize = nothing
    )

    # 1. Check whether the orbitals are Kramer-paired.
    # 2. Check whether the location of all orbitals are at the origin (for periodicity of the Hamiltonian in momentum space)
    # 3. Check whether the the hoppings can be Kramer-paired

    # 4. Construct Hamiltonian as function of momentum
    # 5. Check whether the H(0) and H(G) are the same
    # 6. Check whether H(k) = U ⋅ H(-k) ⋅ U⁻¹
    error("unimplemented")
end



"""
generate k-space grid
(which has (2 n1, 2 n2) points TOTAL in the Brillouin zone)

# Example

When n1 = 4, n2 = 3, this function returns an `OrderedDict` that represents the following structure

```
i2|
  |
5 | -i -i -i -i -i -i -i -i
4 | -i -i -i -i -i -i -i -i
3 | 0h +h +h +h 0h -h -h -h
2 | +i +i +i +i +i +i +i +i
1 | +i +i +i +i +i +i +i +i
0 | 0z +z +z +z 0z -z -z -z
--+----------------------------
  |  0  1  2  3  4  5  6  7  i1
```

where 0z, +z, -z are represented respectively by `:TRIZERO`, `:POSZERO`, and `:NEGZERO`,
and   0h, +h, -h by `:TRIHALF`, `:POSHALF`, and `:NEGHALF`,
and   +i, -i by `:POSINT`, `:NEGINT`.
"""
function timereversalindexgrid(n1 ::Integer, n2 ::Integer)
    pointtypes = OrderedDict{Vector{Int}, Tuple{Symbol, Vector{Int}}}()
    for i2 in 0:(2*n2-1), i1 in 0:(2*n1-1)
        j1 = mod(-i1, 2*n1)
        j2 = mod(-i2, 2*n2)
        if i2 == 0
            if (i1, i2) == (j1, j2)
                pointtypes[[i1,i2]] = (:TRIZERO, [i1,i2])
            elseif haskey(pointtypes, [j1,j2])
                pointtypes[[i1,i2]] = (:NEGZERO, [j1,j2])
            else
                pointtypes[[i1,i2]] = (:POSZERO, [i1,i2])
            end
        elseif i2 == n2
            if (i1, i2) == (j1, j2)
                pointtypes[[i1,i2]] = (:TRIHALF, [i1,i2])
            elseif haskey(pointtypes, [j1,j2])
                pointtypes[[i1,i2]] = (:NEGHALF, [j1,j2])
            else
                pointtypes[[i1,i2]] = (:POSHALF, [i1,i2])
            end
        else
            if (i1, i2) == (j1, j2)
                @assert(false)
            elseif haskey(pointtypes, [j1,j2])
                pointtypes[[i1,i2]] = (:NEGINT, [j1,j2])
            else
                pointtypes[[i1,i2]] = (:POSINT, [i1,i2])
            end
        end
    end
    return pointtypes
end


function degenerategroups(values::AbstractVector{<:Real};
                          tolerance::Real = sqrt(eps(Float64))) ::Vector{Vector{Int}}
    @assert(tolerance > 0)

    n = length(values)

    diffs = values[2:end] - values[1:end-1]

    dgs = Vector{Int}[]
    cg = Int[1]

    for i in 2:n
        if diffs[i-1] < tolerance
            push!(cg, i)
        else
            push!(dgs, cg)
            cg = Int[i]
        end
    end

    if !isempty(cg)
        push!(dgs, cg)
    end
    @debug "Degenerate groups"
    @debug "  values: $values"
    @debug "  groups: $dgs"
    return dgs
end


"""

Get the phase of hamiltonian function.
    Returns nan if no phase

"""
function getnambuphase(
    unitcell ::UnitCell{O},
    Deltafunction ::Function,
    timereversalmatrix ::AbstractMatrix{<:Number},
    n1 ::Integer, n2 ::Integer,
    ) where {O}

    squareunitcell = squarify(unitcell)
    kgrid = momentumgrid(squareunitcell, [2n1, 2n2])
    igrid = timereversalindexgrid(n1, n2)

    Δgrid = Dict()
    for (idx, (t, idx2)) in igrid
        k = kgrid[(i+1 for i in idx)...]
        Δ = Deltafunction(k)
        Δgrid[idx] = Δ
    end

    phasegrid = Dict()
    for (idx, (t, idx2)) in igrid
        Δ1 = Δgrid[idx]
        Δ2 = Δgrid[idx2]
        Δ1c = conj(timereversalmatrix * Δ * timereversalmatrix')
        Δ1cf = [Iterators.flatten(Δ1c)...]
        Δ2f  = [Iterators.flatten(Δ2)...]

        Δ1cn = norm(Δ1cf)
        Δ2n  = norm(Δ2f)

        @assert( isapprox(Δ1cn, Δ2n; atol=tolerance) )

        phase = if Δ1cn > tolerance && Δ2n > tolerance
            phi = dot(Δ1cf, Δ2f)
            phi / abs(phi)
        else
            NaN
        end

        phasegrid[idx] = phi

        #=
        if t == :TRIZERO || t == :TRIHALF

        elseif t == :POSZERO || t == :POSHALF || t == :POSINT

        elseif t == :NEGZERO || t == :NEGHALF

        else
            @error("Unknown type $t")
        end
        =#
    end


end





function kramerpairup!(eigenvalues::AbstractVector{<:Real},
                      eigenvectors::AbstractMatrix{<:Number},
                      theta::AbstractMatrix{<:Number};
                      tolerance::Real = sqrt(eps(Float64)))

    n = length(eigenvalues)
    n1, n2 = size(eigenvectors)
    @assert(n == n1 && n == n2)

    dgs = degenerategroups(eigenvalues; tolerance=tolerance)
    @assert(all(mod(length(x), 2) == 0 for x in dgs), "eigenvalues should come in pairs ($eigenvalues)")

    for dg in dgs
        kramerpairupblock!(view(eigenvectors, :, dg), theta; tolerance=tolerance)
    end


end


function kramerpairupblock!(a ::AbstractMatrix{C1},
                            theta ::AbstractMatrix{C2};
                            tolerance::Real = sqrt(eps(Float64))) where {C1<:Number, C2<:Number}
    @assert(promote_type(C1, C2) == C1)

    n, m = size(a)
    @debug "n = $n, m = $m"

    vectorqueue = Vector{C1}[a[:,i] for i in 1:m]
    vectorout = Vector{C1}[]
    @assert(mod(m, 2) == 0)


    for i in 1:(div(m, 2))
        sort!(vectorqueue, by=norm, rev=true)

        largestnorm = norm(vectorqueue[1])
        for (j, vec) in enumerate(vectorqueue)
          vectorqueue[j] /= largestnorm
        end

        psi1 = vectorqueue[1]
        psi2 = theta * conj(psi1)
        @assert(isapprox(dot(psi1, psi2), 0; atol=tolerance), "psi1 not orthogonal to psi2")
        normalize!(psi1)
        normalize!(psi2)

        push!(vectorout, psi1)
        push!(vectorout, psi2)

        newvectornorm = norm(vectorqueue[1])
        newvectorqueue = Vector{C1}[]
        for j in 2:length(vectorqueue)
            vec = vectorqueue[j]
            for psi in vectorout
                vec -= psi * dot(psi, vec)
            end
            #newvec = vec - psi1 * dot(psi1, vec) - psi2 * dot(psi2, vec)
            push!(newvectorqueue, vec)
        end
        vectorqueue = newvectorqueue
    end

    let
        norms = [norm(vec) for vec in vectorqueue]
        if maximum(norms) > tolerance
            @warn "norms of remaining vectors exceeds tolerance"
            @warn "  norms = $(norms)"
            @warn "  tolerance = $(tolerance)"
        end
    end
    #@assert(all(isapprox(norm(vec), 0; atol=tolerance) for vec in vectorqueue), "Remaining \"vectors\" in vectorqueue should span zero dimension")
    for i in 1:m
        a[:, i] = vectorout[i]
    end
end
