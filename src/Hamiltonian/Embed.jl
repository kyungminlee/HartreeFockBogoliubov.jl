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
struct HoppingDiagonal{R<:Real}
  amplitude ::R
  i ::Int64
  ri ::CarteCoord
end

"""
    HoppingDiagonal
"""
function HoppingDiagonal(uc ::UnitCell{O},
                         hopspec ::Spec.HoppingDiagonal{R}) where {O, R<:Real}
  i = hopspec.i
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  return HoppingDiagonal{R}(hopspec.amplitude, i, ri)
end


"""
    HoppingOffdiagonal
"""
struct HoppingOffdiagonal{C<:Number}
  amplitude ::C
  i ::Int64
  j ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
end


"""
    HoppingOffdiagonal
"""
function HoppingOffdiagonal(uc ::UnitCell{O},
                            hopspec ::Spec.HoppingOffdiagonal{C}) where {O, C<:Number}
  (i, j) = (hopspec.i, hopspec.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)

  v = hopspec.amplitude
  @assert(i <= j, "ordering of i and j should have been taken care of by Spec")
  return HoppingOffdiagonal{C}(v, i, j, ri, rj)
end


"""
    InteractionDiagonal
"""
struct InteractionDiagonal{R<:Real}
  amplitude ::R
  i ::Int64
  j ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
end


"""
    InteractionDiagonal
"""
function InteractionDiagonal(uc::UnitCell,
                             hopspec::Spec.InteractionDiagonal{R}) where {R<:Real}
  (i, j) = (hopspec.i, hopspec.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)
  v = hopspec.amplitude
  @assert(i <= j, "ordering of i and j should have been taken care of by Spec")
  return InteractionDiagonal{R}(v, i, j, ri, rj)
end


"""
    InteractionOffdiagonal
"""
struct InteractionOffdiagonal{C<:Number}
  amplitude ::C
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
function InteractionOffdiagonal(uc::UnitCell{O},
                                hopspec::Spec.InteractionOffdiagonal{C}) where {O, C<:Number}
  (i, j) = (hopspec.i, hopspec.j)
  (k, l) = (hopspec.k, hopspec.l)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)
  rk = fract2carte(uc, getorbitalcoord(uc, hopspec.k) + hopspec.Rk)
  rl = fract2carte(uc, getorbitalcoord(uc, hopspec.l) + hopspec.Rl)
  v = hopspec.amplitude
  @assert(i <= j && k <= l && i <= k)
  return InteractionOffdiagonal{C}(v, i, j, k, l, ri, rj, rk, rl)
end


"""
    embed
"""
function embed(uc ::UnitCell{O},
               hopspec ::Spec.HoppingDiagonal{R}) where {O, R<:Real}
  return HoppingDiagonal(uc, hopspec)
end


"""
    embed
"""
function embed(uc ::UnitCell{O},
               hopspec ::Spec.HoppingOffdiagonal{C}) where {O, C<:Number}
  return HoppingOffdiagonal(uc, hopspec)
end


"""
    embed
"""
function embed(uc::UnitCell{O},
               hopspec::Spec.InteractionDiagonal{R}) where {O, R<:Real}
  return InteractionDiagonal(uc, hopspec)
end


"""
    embed
"""
function embed(uc::UnitCell{O},
               hopspec::Spec.InteractionOffdiagonal{C}) where {O, C<:Number}
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
mutable struct Hamiltonian{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
end


"""
    Hamiltonian
"""
function Hamiltonian(spec ::Spec.Hamiltonian{O}) where {O}
  unitcell = spec.unitcell
  hoppings = [embed(unitcell, hop) for hop in spec.hoppings]
  interactions = [embed(unitcell, inter) for inter in spec.interactions]
  Hamiltonian{O}(unitcell, hoppings, interactions)
end

end # Embed
