module Generator

using ..Lattice

import ..Spec
import ..Embed


function generatefast(uc ::UnitCell, hop ::Embed.HoppingDiagonal)
  ndim = dimension(uc)
  norb = numorbital(uc)
  v = hop.amplitude
  i = hop.i
  return (momentum ::Vector{Float64}, out::Array{Complex128, 2}) -> begin
    @assert(size(momentum) == (ndim,))
    @assert(size(out) == (norb, norb))
    out[i,i] += v
    return out
  end
end


function generatefast(uc ::UnitCell, hop::Embed.HoppingOffdiagonal)
  ndim = dimension(uc)
  norb = numorbital(uc)
  v = hop.amplitude
  (i, j) = (hop.i, hop.j)
  (ri, rj) = (hop.ri, hop.rj)
  rji = fract2carte(uc, rj) - fract2carte(uc, ri)

  return (momentum ::Vector{Float64}, out::Array{Complex128, 2}) -> begin
    @assert(size(momentum) == (ndim,))
    @assert(size(out) == (norb, norb))
    phase = exp(1im * dot(momentum, rji))
    out[i,j] += v * phase
    out[j,i] += conj(v * phase)
    return out
  end
end


function generatefast(uc ::UnitCell, hops ::Vector{Embed.Hopping})
  ndim = dimension(uc)
  norb = numorbital(uc)
  funcs = [generatefast(uc, hop) for hop in hops]
  return (momentum ::Vector{Float64}, out::Array{Complex128, 2}) -> begin
    for func in funcs
      func(momentum, out)
    end
    return out
  end
end


function generatehoppingfast(hamiltonian ::Embed.Hamiltonian)
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end


function generatefast(uc ::UnitCell, hopspec::Spec.HoppingDiagonal)
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


function generatefast(uc ::UnitCell, hopspec::Spec.HoppingOffdiagonal)
  hopembed = Embed.embed(uc, hopspec)
  return generatefast(uc, hopembed)
end


function generatefast(uc ::UnitCell, hopspecs ::Vector{Spec.Hopping})
  hopembeds = Embed.Hopping[Embed.embed(uc, hopspec) for hopspec in hopspecs]
  return generatefast(uc, hopembeds)
end


function generatehoppingfast(hamiltonian ::Spec.Hamiltonian)
  return generatefast(hamiltonian.unitcell, hamiltonian.hoppings)
end

end
