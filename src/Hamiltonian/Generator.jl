module Generator

using ..Lattice

import ..Spec
import ..Embed

"""
    generatefast
"""
function generatefast(uc ::UnitCell{O},
                      hop ::Embed.HoppingDiagonal{R}) where {O, R<:Real}
  ndim = dimension(uc)
  norb = numorbital(uc)
  v = hop.amplitude
  i = hop.i
  return (momentum ::AbstractVector{Float64}, out::AbstractArray{Complex128, 2}) -> begin
    @assert(size(momentum) == (ndim,))
    @assert(size(out) == (norb, norb))
    out[i,i] += v
    return out
  end
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O},
                      hop::Embed.HoppingOffdiagonal{C}) where {O, C<:Number}
  ndim = dimension(uc)
  norb = numorbital(uc)
  v = hop.amplitude
  (i, j) = (hop.i, hop.j)
  (ri, rj) = (hop.ri, hop.rj)
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
                      hops ::AbstractVector{Embed.Hopping}) where {O}
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
    generatehoppingfast
"""
function generatehoppingfast(hamiltonian ::Embed.Hamiltonian{O}) where {O}
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O}, hopspec::Spec.HoppingDiagonal) where {O}
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O}, hopspec::Spec.HoppingOffdiagonal) where {O}
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O},
                      hopspecs ::AbstractVector{Spec.Hopping}) where {O}
  hopembeds = Embed.Hopping[Embed.embed(uc, hopspec) for hopspec in hopspecs]
  return generatefast(uc, hopembeds)
end


"""
    generatefast
"""
function generatehoppingfast(hamiltonian ::Spec.Hamiltonian{O}) where {O}
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end

end
