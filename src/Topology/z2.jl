import DataStructures: OrderedDict

function iscanonicaltimereversalinvariant(
    unitcell::UnitCell,
    timereversalmatrix::AbstractMatrix{<:Number},
    hoppings::AbstractVector{Hopping};
    tol::Real = sqrt(eps(Float64)),
    momentumgridsize
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
    pointtypes = OrderedDict{Vector{Int64}, Tuple{Symbol, Vector{Int64}}}()
    for i2 in 0:(2*n2-1), i1 in 0:(2*n1 - 1)
        j1 = mod(-i1, 2*n1)
        j2 = mod(-i2, 2*n2)
        if i2 == 0 || i2 == n2
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

#=
function timereversalindexgrid(shape::AbstractVector{<:Integer})
ranges = [0:(2*n-1) for n in shape]
hypercubicgrid = map((x) -> [x...], Base.product(ranges...))

pointtypes = OrderedDict{Vector{Int64}, Tuple{Symbol, Vector{Int64}}}()

for i in hypercubicgrid
j = [mod(-x, 2*n) for (x, n) in zip(i, shape)]
if (i[end] == 0 || i[end] == shape[end])
if i == j
pointtypes[i] = (:TRI, i)
elseif haskey(pointtypes, j)
pointtypes[i] = (:NEG, j)
else
pointtypes[i] = (:POS, i)
end
else
if i == j
@assert(False)
elseif haskey(pointtypes, j)
pointtypes[i] = (:NEGINT, j)
else
pointtypes[i] = (:POSINT, i)
end
end
end
return pointtypes
end
=#

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
The Z2 index
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
            all(isapprox(hk0, hk1)) && all(isapprox(hk0, hk2))
        end,
        "Hamiltonian at reciprocal lattice vectors should be the same as at 0.")

    hk = zeros(Complex128, norb, norb)
    for (idx, (t, idx2)) in igrid
        if t == :TRIZERO || t == :TRIHALF
            k = kgrid[(i+1 for i in idx)...]
            hk[:] = 0.0
            hkgen(k, hk)
            u, v = eig(Hermitian(0.5 * (hk + hk')))
            for idxpair in selectpairs
                @assert(isapprox(u[idxpair*2-1], u[idxpair*2]; rtol=rtol, atol=atol))
                v[:, idxpair*2] = timereversal * conj(v[:, idxpair*2-1])    # T = U_T K
            end
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
            t == :TRI || t == :POS || t == :NEG || t == :POSINT
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
            t == :TRI || t == :POS || t == :NEG
        end)

    ψ1 = eigenvectorgrid[[i1 ,0]]
    ψ2 = eigenvectorgrid[[i1p,0]]
    A += angle(det(ψ1' * ψ2))
    end

    for i1 in 0:(2*n1-1)
        @assert(let
            (t,_) = igrid[[i1,n2]]
            t == :TRI || t == :POS || t == :NEG
        end)
        i1p = mod(i1 + 1, 2*n1)
        ψ1 = eigenvectorgrid[[i1 ,n2]]
        ψ2 = eigenvectorgrid[[i1p,n2]]
        A -= angle(det(ψ1' * ψ2))
    end
    z2indexreal = (F-A) / (2π)
    z2indexint = round(z2indexreal)
    @assert(abs(z2indexreal - z2indexint) < atol)
    return mod(z2indexint, 2)
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
