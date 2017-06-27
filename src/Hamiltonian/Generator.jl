module Generator

using ..Lattice

using ..Spec


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O}, hop::SpecHoppingDiagonal) where {O}
  ndim = dimension(uc)
  norb = numorbital(uc)
  v = hop.amplitude
  i = hop.i

  function(momentum ::AbstractVector{Float64}, out::AbstractArray{Complex128, 2})
    @assert(size(momentum) == (ndim,))
    @assert(size(out) == (norb, norb))
    out[i,i] += v
    return out
  end
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O}, hop::SpecHoppingOffdiagonal) where {O}
  ndim = dimension(uc)
  norb = numorbital(uc)
  v = hop.amplitude
  (i, j) = (hop.i, hop.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hop.i) + hop.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hop.j) + hop.Rj)
  rij = rj - ri

  function(momentum ::AbstractVector{Float64}, out::AbstractArray{Complex128, 2})
    @assert(size(momentum) == (ndim,))
    @assert(size(out) == (norb, norb))
    phase = exp(1im * dot(momentum, rij))
    out[i,j] += v * phase
    out[j,i] += conj(v * phase)
    return out
  end
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O},
                      hops ::AbstractVector{SpecHopping}) where {O}
  ndim = dimension(uc)
  norb = numorbital(uc)
  funcs = [generatefast(uc, hop) for hop in hops]

  function(momentum ::AbstractVector{Float64}, out::AbstractArray{Complex128, 2})
    for func in funcs
      func(momentum, out)
    end
    return out
  end
end


"""
    generatefast
"""
function generatehoppingfast(hamiltonian ::SpecHamiltonian{O}) where {O}
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end

end
