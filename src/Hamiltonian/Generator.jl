module Generator

using ..Lattice

using ..Spec
using ..Embed

"""
    generatefast
"""
function generatefast(uc ::UnitCell{O},
                      hop ::EmbedHoppingDiagonal{R}) where {O, R<:Real}
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
function generatefast(uc ::UnitCell{O},
                      hop::EmbedHoppingOffdiagonal{C}) where {O, C<:Number}
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
                      hops ::AbstractVector{EmbedHopping}) where {O}
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
function generatehoppingfast(hamiltonian ::EmbedHamiltonian{O}) where {O}
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O}, hopspec::SpecHoppingDiagonal) where {O}
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O}, hopspec::SpecHoppingOffdiagonal) where {O}
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


"""
    generatefast
"""
function generatefast(uc ::UnitCell{O},
                      hopspecs ::AbstractVector{SpecHopping}) where {O}
  hopembeds = Embed.EmbedHopping[Embed.embed(uc, hopspec) for hopspec in hopspecs]
  return generatefast(uc, hopembeds)
end


"""
    generatefast
"""
function generatehoppingfast(hamiltonian ::SpecHamiltonian{O}) where {O}
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end

end
