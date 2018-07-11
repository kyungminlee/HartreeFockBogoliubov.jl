if VERSION < v"0.7-"
    using MicroLogging
end
import DataStructures: OrderedDict

"""
z2index

Compute Z2 index of time-reversal-invariant Hamiltonian.

# Arguments
* `uc::UnitCell{O}`
* `hops::AbstractVector{Hopping}`
* `timereversal::AbstractMatrix`
* `n1 ::Integer`
* `n2 ::Integer`
* `selectpairs::AbstractVector{<:Integer}`

# Optional Arguments
* `tol ::Real = sqrt(eps(Float64))`

# Returns
(The Z2 index,  max| Hₖ - T⁻¹HₖT | for k in TRIMs)
"""
function z2index(uc::UnitCell{O},
                 hops::AbstractVector{Hopping},
                 timereversal::AbstractMatrix,
                 n1 ::Integer,
                 n2 ::Integer,
                 selectpairs::AbstractVector{<:Integer},
                 ;
                 rtol ::Real = sqrt(eps(Float64)),
                 atol ::Real = sqrt(eps(Float64)),
                 ) where {O}

    squareuc = squarify(uc)
    norb = numorbital(squareuc)
    @assert(mod(norb, 2) == 0)
    hkgen = Generator.generatefast(squareuc, hops)

    kgrid = momentumgrid(squareuc, [n1*2, n2*2])
    igrid = timereversalindexgrid(n1, n2)

    eigenvaluegrid = Dict()
    eigenvectorgrid = Dict()

    selectbands = Int[]
    for idxpair in selectpairs
        @assert(1 <= idxpair <= norb ÷ 2)
        push!(selectbands, idxpair*2-1)
        push!(selectbands, idxpair*2)
    end

    # check time reversal
    @assert(all(isapprox.(timereversal + timereversal.', 0.0; rtol=rtol, atol=atol)),
            "timereversalmatrix need to be antisymmetric")
    @assert(isapprox(timereversal*timereversal', eye(norb); rtol=rtol, atol=atol),
            "timereversalmatrix need to be unitary")

    @assert(let
            hk0 = zeros(Complex128, norb, norb)
            hk1 = zeros(Complex128, norb, norb)
            hk2 = zeros(Complex128, norb, norb)

            hkgen([0.0, 0.0], hk0)
            hkgen(squareuc.reciprocallatticevectors[:,1], hk1)
            hkgen(squareuc.reciprocallatticevectors[:,2], hk2)

            all(isapprox(hk0, hk1; rtol=rtol, atol=atol)) && all(isapprox(hk0, hk2; rtol=rtol, atol=atol))
        end,
        "Hamiltonian at reciprocal lattice vectors should be the same as that at 0.")

    maxDiffTimeReversal = 0.0

    hk = zeros(Complex128, norb, norb)

    # Check Time Reversal First
    for (idx, (t, idx2)) in igrid
        if t == :TRIZERO || t == :TRIHALF
            k = kgrid[(i+1 for i in idx)...]
            hk[:] = 0.0
            hkgen(k, hk)
            hk2 = timereversal' * hk * timereversal

            maxdiff = maximum(abs.((hk2).' - hk))
            maxDiffTimeReversal = max(maxDiffTimeReversal, maxdiff)
        end
    end

    if maxDiffTimeReversal > atol
        @warn "Hamiltonian is not time reversal symmetric ( max| Hₖ - T⁻¹HₖT | = $maxDiffTimeReversal )"
        return ( NaN, maxDiffTimeReversal )
    end

    for (idx, (t, idx2)) in igrid
        if t == :TRIZERO || t == :TRIHALF
            k = kgrid[(i+1 for i in idx)...]
            hk[:] = 0.0
            hkgen(k, hk)
            u, v = eig(Hermitian(0.5 * (hk + hk')))
            @debug("Eigenvalues at TRIM ($k) = $u")

            kramerpairup!(u, v, timereversal; tolerance=max(atol, rtol))
            eigenvaluegrid[idx] = u[selectbands]
            eigenvectorgrid[idx] = v[:, selectbands]
        elseif t == :POSZERO || t == :POSHALF || t == :POSINT
            k = kgrid[(i+1 for i in idx)...]
            hk[:] = 0.0
            hkgen(k, hk)
            u, v = eig(Hermitian(0.5 * (hk + hk')))
            eigenvaluegrid[idx] = u[selectbands]
            eigenvectorgrid[idx] = v[:, selectbands]
        end
    end

    for (idx, (t, idx2)) in igrid
        if t == :NEGZERO || t == :NEGHALF
            eigenvalues = eigenvaluegrid[idx2]
            eigenvectors = eigenvectorgrid[idx2]
            eigenvaluegrid[idx] = eigenvalues
            eigenvectorgrid[idx] = timereversal * conj(eigenvectors)
        end
    end

    # Eigenvectorgrid done. Now compute F and A

    F = 0.0
    for i2 in 0:(n2-1), i1 in 0:(2*n1-1)
        @assert(let
            (t,_) = igrid[[i1,i2]]
            t == :TRIZERO || t == :POSZERO || t == :NEGZERO || t == :POSINT
        end)
        i1p = mod(i1 + 1, 2*n1)
        i2p = i2 + 1
        ψ1 = eigenvectorgrid[[i1 ,i2 ]]
        ψ2 = eigenvectorgrid[[i1p,i2 ]]
        ψ3 = eigenvectorgrid[[i1p,i2p]]
        ψ4 = eigenvectorgrid[[i1 ,i2p]]
        F += angle(det(ψ1' * ψ2) * det(ψ2' * ψ3) * det(ψ3' * ψ4) * det(ψ4' * ψ1))
    end

    A = 0.0
    for i1 in 0:(2*n1-1)
        i1p = mod(i1 + 1, 2*n1)

        @assert(let
                (t,_) = igrid[[i1,0]]
                t == :TRIZERO || t == :POSZERO || t == :NEGZERO
            end,
            "$t has the wrong type"
            )

    ψ1 = eigenvectorgrid[[i1 ,0]]
    ψ2 = eigenvectorgrid[[i1p,0]]
    A += angle(det(ψ1' * ψ2))
    end

    for i1 in 0:(2*n1-1)
        @assert(let
                (t,_) = igrid[[i1,n2]]
                t == :TRIHALF || t == :POSHALF || t == :NEGHALF
            end,
            "$t has the wrong type"
            )
        i1p = mod(i1 + 1, 2*n1)
        ψ1 = eigenvectorgrid[[i1 ,n2]]
        ψ2 = eigenvectorgrid[[i1p,n2]]
        A -= angle(det(ψ1' * ψ2))
    end
    z2indexreal = (F-A) / (2π)
    z2indexint = round(Int, z2indexreal)
    if ! ( abs(z2indexreal - z2indexint) < atol )
        @warn "$z2indexreal is not close to integer"
    end

    @info "A/2π = $(A/2π)"
    @info "F/2π = $(F/2π)"


    return ( mod(round(z2indexreal, precision(z2indexreal) ÷ 2, 2), 2), maxDiffTimeReversal )
    #return mod(z2indexint, 2)
end

#=
function z2invariant(uc::UnitCell{O},
hops::AbstractVector{Hopping},
n1::Integer,
n2::Integer,
selectband::AbstractVector{<:Integer}) where {O}
@assert(n1 > 0 && n1 % 2 == 0, "n1 should be a positive even number")
@assert(n2 > 0 && n2 % 2 == 0, "n2 should be a positive even number")

@assert(length(selectband) % 2 == 0, "selectband should account for time reversal")

squareuc = squarify(uc)
igrid = triindexgrid(n1, n2)
eg

end
=#
