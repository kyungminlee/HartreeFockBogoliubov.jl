module Generator

using ..Lattice

import ..Spec
import ..Embed

"""
    generatefast
"""
function generatefast{T}(uc ::UnitCell{T}, hop ::Embed.HoppingDiagonal)
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
function generatefast{T}(uc ::UnitCell{T}, hop::Embed.HoppingOffdiagonal)
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
function generatefast{T}(uc ::UnitCell{T}, hops ::Vector{Embed.Hopping})
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
function generatehoppingfast{T}(hamiltonian ::Embed.Hamiltonian{T})
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end


"""
    generatefast
"""
function generatefast{T}(uc ::UnitCell{T}, hopspec::Spec.HoppingDiagonal)
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


"""
    generatefast
"""
function generatefast{T}(uc ::UnitCell{T}, hopspec::Spec.HoppingOffdiagonal)
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


"""
    generatefast
"""
function generatefast{T}(uc ::UnitCell{T}, hopspecs ::AbstractVector{Spec.Hopping})
  hopembeds = Embed.Hopping[Embed.embed(uc, hopspec) for hopspec in hopspecs]
  return generatefast(uc, hopembeds)
end


"""
    generatefast
"""
function generatehoppingfast{T}(hamiltonian ::Spec.Hamiltonian{T})
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end

end
