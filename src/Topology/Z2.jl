export squarify

"""
generate k-space grid
(which has (2 n1, 2 n2) points TOTAL in the Brillouin zone)
"""
function triindexgrid(n1 ::Integer, n2 ::Integer)
    pointtypes = OrderedDict{Vector{Int64}, Tuple{Symbol, Vector{Int64}}}()
    for i2 in 0:(2*n2-1), i1 in 0:(2*n1 - 1)
        j1 = mod(-i1, 2*n1)
        j2 = mod(-i2, 2*n2)
        if (i1, i2) == (j1, j2)
            pointtypes[[i1,i2]] = (:TRI, [i1,i2])
        else
            if haskey(pointtypes, [j1,j2])
                pointtypes[[i1,i2]] = (:NEG, [j1,j2])
            else
                pointtypes[[i1,i2]] = (:POS, [i1,i2])
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
        if i == j
            pointtypes[i] = (:TRI, i)
        elseif haskey(pointtypes, j)
            pointtypes[i] = (:NEG, j)
        else
            pointtypes[i] = (:POS, i)
        end
    end
    return pointtypes
end

function eigengrid(func::Function,
                     timereversal::Function,
                     reciprocallatticevectors::Vector{Float64},
                     n1 ::Integer,
                     n2 ::Integer,
                     )
  igrid = triindexgrid(n1, n2)
  eg = Dict()

  for (idx, (t, idx2)) in igrid
    if t == :TRI || t == :POS
      k = reciprocallatticevectors * [idx[1] / n1, idx[2] / n2]
      hk = func(k)
      eg[k1] = eig(Hermitian(0.5 * (hk + hk')))
    end
  end

  for (idx, (t, idx2)) in igrid
    if t == :NEG
      (eigenvalues, eigenvectors) = eg[idx2]
      eg[k1] = (eigenvalues, timereversal(eigenvectors))
    end
  end

  return eg
end


"""
    squarify

  # Arguments
  * `uc::Lattice.UnitCell{O}`
"""
function squarify(uc::Lattice.UnitCell{O}) where {O}
  newdimension = dimension(uc)
  newlatticevectors = eye(newdimension)

  origin = FractCoord(zeros(Int64, newdimension), zeros(Float64, newdimension))
  
  newuc = newunitcell(newlatticevectors; OrbitalType=O)
  for (orbname, fc) in uc.orbitals
    addorbital!(newuc, orbname, origin)
  end
  return newuc
end

"""
    squarify

  # Arguments
  * `uc::Spec.FullHamiltonian{O}`
"""
function squarify(ham::Spec.FullHamiltonian{O}) where {O}
  newuc = squarify(ham.unitcell)
  return FullHamiltonian{O}(newuc, ham.hoppings, ham.interactions)
end

"""
    squarify

  # Arguments
  * `uc::HFB.HFBHamiltonian{O}`
"""
function squarify(ham::HFB.HFBHamiltonian{O}) where {O}
  newuc = squarify(ham.unitcell)
  return HFBHamiltonian{O}(newuc,
                           ham.hoppings,
                           ham.particle_hole_interactions,
                           ham.particle_particle_interactions)
end