"""
    Submodule `Embed`

From `Spec` for Hamiltonian in periodic system, create an `Embed`


"""
module Embed

using ..Lattice
import ..Spec


"""
    HoppingDiagonal
"""
immutable HoppingDiagonal
  amplitude ::Real
  i ::Int64
  ri ::CarteCoord
end


"""
    HoppingDiagonal
"""
function HoppingDiagonal(uc ::UnitCell, hopspec ::Spec.HoppingDiagonal)
  i = hopspec.i
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  return HoppingDiagonal(hopspec.amplitude, i, ri)
end


"""
    HoppingOffdiagonal
"""
immutable HoppingOffdiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
end


"""
    HoppingOffdiagonal
"""
function HoppingOffdiagonal(uc ::UnitCell, hopspec ::Spec.HoppingOffdiagonal)
  (i, j) = (hopspec.i, hopspec.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)

  v = hopspec.amplitude
  @assert(i <= j, "ordering of i and j should have been taken care of by Spec")
  return HoppingOffdiagonal(v, i, j, ri, rj)
end


"""
    InteractionDiagonal
"""
immutable InteractionDiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
end


"""
    InteractionDiagonal
"""
function InteractionDiagonal(uc ::UnitCell, hopspec ::Spec.InteractionDiagonal)
  (i, j) = (hopspec.i, hopspec.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)
  v = hopspec.amplitude
  @assert(i <= j, "ordering of i and j should have been taken care of by Spec")
  return InteractionDiagonal(v, i, j, ri, rj)
end


"""
    InteractionOffdiagonal
"""
immutable InteractionOffdiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  k ::Int64
  l ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
  rk ::CarteCoord
  rl ::CarteCoord
end


"""
    InteractionOffdiagonal
"""
function InteractionOffdiagonal(uc ::UnitCell, hopspec ::Spec.InteractionOffdiagonal)
  (i, j) = (hopspec.i, hopspec.j)
  (k, l) = (hopspec.k, hopspec.l)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)
  rk = fract2carte(uc, getorbitalcoord(uc, hopspec.k) + hopspec.Rk)
  rl = fract2carte(uc, getorbitalcoord(uc, hopspec.l) + hopspec.Rl)
  v = hopspec.amplitude
  @assert(i <= j && k <= l && i <= k)
  return InteractionOffdiagonal(v, i, j, k, l, ri, rj, rk, rl)
end


"""
    embed
"""
function embed{T}(uc ::UnitCell{T}, hopspec ::Spec.HoppingDiagonal)
  return HoppingDiagonal(uc, hopspec)
end


"""
    embed
"""
function embed{T}(uc ::UnitCell{T}, hopspec ::Spec.HoppingOffdiagonal)
  return HoppingOffdiagonal(uc, hopspec)
end


"""
    embed
"""
function embed{T}(uc ::UnitCell{T}, hopspec ::Spec.InteractionDiagonal)
  return InteractionDiagonal(uc, hopspec)
end


"""
    embed
"""
function embed{T}(uc ::UnitCell{T}, hopspec ::Spec.InteractionOffdiagonal)
  return InteractionOffdiagonal(uc, hopspec)
end


"""
    Hopping
"""
const Hopping = Union{HoppingDiagonal, HoppingOffdiagonal}


"""
    Interaction
"""
const Interaction = Union{InteractionDiagonal, InteractionOffdiagonal}


"""
    Hamiltonian
"""
type Hamiltonian{T}
  unitcell ::UnitCell{T}
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
end


"""
    Hamiltonian
"""
function Hamiltonian{T}(spec ::Spec.Hamiltonian{T})
  unitcell = spec.unitcell
  hoppings = [embed(unitcell, hop) for hop in spec.hoppings]
  interactions = [embed(unitcell, inter) for inter in spec.interactions]
  Hamiltonian(unitcell, hoppings, interactions)
end

end # Embed
