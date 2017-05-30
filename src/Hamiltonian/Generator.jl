module Generator

using ..Lattice

import ..Spec
import ..Embed

function generate(uc ::UnitCell, hops ::Vector{Embed.Hopping})
  ndim = dimension(uc)
  norb = numorbital(uc)
  funcs = [generate(uc, hop) for hop in hops]
  return (momentum ::Vector{Float64}, out::Array{Complex128, 2}) -> begin
    for func in funcs
      func(momentum, out)
    end
    return out
  end
end

function generate(uc ::UnitCell, hop ::Embed.HoppingDiagonal)
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

function generate(uc ::UnitCell, hop::Embed.HoppingOffdiagonal)
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
    #@show i, j, phase
    out[i,j] += v * phase
    out[j,i] += conj(v * phase)
    return out
  end
end




end