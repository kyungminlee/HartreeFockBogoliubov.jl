export triindexgrid

import DataStructures: OrderedDict

"""
generate k-space grid
(which has (2 n1, 2 n2) points TOTAL in the Brillouin zone)
"""
function triindexgrid(n1 ::Integer, n2 ::Integer)
    pointtypes = OrderedDict{Vector{Int64}, Tuple{Symbol, Vector{Int64}}}()
    for i2 in 0:(2*n2-1), i1 in 0:(2*n1 - 1)
        j1 = mod(-i1, 2*n1)
        j2 = mod(-i2, 2*n2)
        if i2 == 0 || i2 == n2
            if (i1, i2) == (j1, j2)
                pointtypes[[i1,i2]] = (:TRI, [i1,i2])
            elseif haskey(pointtypes, [j1,j2])
                pointtypes[[i1,i2]] = (:NEG, [j1,j2])
            else
                pointtypes[[i1,i2]] = (:POS, [i1,i2])
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

function triindexgrid(shape::AbstractVector{<:Integer})
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

function z2invariant(uc::UnitCell{O},
                     hops::AbstractVector{Hopping},
                     timereversal::AbstractMatrix,
                     n1 ::Integer,
                     n2 ::Integer,
                     selectpairs::AbstractVector{<:Integer},
                     ;
                     rtol ::Real = 0,
                     atol ::Real = sqrt(eps(Float64)),
                     ) where {O}


  squareuc = squarify(uc)
  igrid = triindexgrid(n1, n2)
  eigenvaluegrid = Dict()
  eigenvectorgrid = Dict()

  selectbands = Int[]
  for idxpair in selectpairs
    push!(selectbands, idxpair*2-1)
    push!(selectbands, idxpair*2)
  end

  #@show timereversal
  #@show igrid
  #@show selectpairs
  #@show selectbands

  norb = numorbital(squareuc)
  hkgen = Generator.generatefast(squareuc, hops)
  kgrid = momentumgrid(squareuc, [n1*2, n2*2])

  # check time reversal
  @assert(all(isapprox.(timereversal + timereversal.', 0.0; rtol=rtol, atol=atol)))
  @assert(isapprox(timereversal*timereversal', eye(norb); rtol=rtol, atol=atol))

  @assert(let
      hk0 = zeros(Complex128, norb, norb)
      hk1 = zeros(Complex128, norb, norb)
      hk2 = zeros(Complex128, norb, norb)
      
      hkgen([0.0, 0.0], hk0)
      hkgen(squareuc.reciprocallatticevectors[:,1], hk1)
      hkgen(squareuc.reciprocallatticevectors[:,2], hk2)
      all(isapprox(hk0, hk1)) && all(isapprox(hk0, hk2))
  end)

  #@show kgrid
  hk = zeros(Complex128, norb, norb)
  for (idx, (t, idx2)) in igrid
    if t == :TRI
      k = kgrid[(i+1 for i in idx)...]
      #@show idx, k, t
      # k = squareuc.reciprocallatticevectors * [idx[1] / (2*n1), idx[2] / (2*n2)]
      hk[:] = 0.0
      hkgen(k, hk)
      u, v = eig(Hermitian(0.5 * (hk + hk')))
      for idxpair in selectpairs
        @assert(isapprox(u[idxpair*2-1], u[idxpair*2]; rtol=rtol, atol=atol))
        #=
        let
            ψ1 = timereversal * conj(v[:, idxpair*2-1]) 
            ψ2 = v[:, idxpair*2]
            (_, idxnonzero) = findmax(abs.(ψ1))
            ψ1 = ψ1 .* conj( ψ1[idxnonzero] / abs(ψ1[idxnonzero]) )
            ψ2 = ψ2 .* conj( ψ2[idxnonzero] / abs(ψ2[idxnonzero]) )
            @assert(isapprox(ψ1, ψ2))
            # NOTE: may not be satisfied when there is additional degeneracy
        end
        =#
        v[:, idxpair*2] = timereversal * conj(v[:, idxpair*2-1])    # T = U_T K
      end
      eigenvaluegrid[idx] = u[selectbands]
      eigenvectorgrid[idx] = v[:, selectbands]
    elseif t == :POS || t == :POSINT
      k = kgrid[(i+1 for i in idx)...]
      #@show idx, k, t
      hk[:] = 0.0
      hkgen(k, hk)
      u, v = eig(Hermitian(0.5 * (hk + hk')))
      eigenvaluegrid[idx] = u[selectbands]
      eigenvectorgrid[idx] = v[:, selectbands]
    end
  end

  for (idx, (t, idx2)) in igrid
    if t == :NEG
      #@show idx, t, idx2
      eigenvalues = eigenvaluegrid[idx2]
      eigenvectors = eigenvectorgrid[idx2]
      eigenvaluegrid[idx] = eigenvalues
      eigenvectorgrid[idx] = timereversal * conj(eigenvectors)
    elseif t == :NEGINT
      #@show idx, t
    end
  end

  #=
  println("Eigenvectorgrid")
  for (k, v) in eigenvectorgrid
    @show k
    display(v)
    println()
  end
  =#

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
    #@show (i1, i2), angle(det(ψ1' * ψ2) * det(ψ2' * ψ3) * det(ψ3' * ψ4) * det(ψ4' * ψ1))
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

    let
        (t,_) = igrid[[i1,0]]
        (tp,_) = igrid[[i1p,0]]
        #@show i1, 0, t, i1p, 0, tp
        #@show ψ1
        #@show ψ2
        #@show '+', angle(det(ψ1' * ψ2))
        #println()
    end
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

    let
        (t,_) = igrid[[i1,0]]
        (tp,_) = igrid[[i1p,0]]
        #@show i1, 0, t, i1p, 0, tp
        #@show ψ1
        #@show ψ2  
        #@show '-', angle(det(ψ1' * ψ2))
        #println()
    end
  end
  #@show mod(F, 4π)
  #@show mod(A, 4π)
  z2indexreal = (F-A) / (2π)
  z2indexint = round(z2indexreal)
  @assert(abs(z2indexreal - z2indexint) < atol)
  #return z2indexint, z2indexreal
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