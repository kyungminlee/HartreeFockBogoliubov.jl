"""
    Submodule `Embed`

From `Spec` for Hamiltonian in periodic system, create an `Embed`


"""

module Embed

using ..Lattice
using ..Spec

export EmbedHoppingDiagonal,
       EmbedHoppingOffdiagonal,
       EmbedInteractionDiagonal,
       EmbedInteractionOffdiagonal
export EmbedHopping, EmbedInteraction
export EmbedHamiltonian
export embed


"""
    HoppingDiagonal

```math
  t c_{i}^{*} c_{i}
```

  # Members
  * `amplitude::R`
  * `i::Int64`: index of orbital i
  * `ri::CarteCoord`: real-space coordinates of orbital i.
"""
struct EmbedHoppingDiagonal{R<:Real}
  amplitude ::R
  i ::Int64
  ri ::CarteCoord
end


"""
    HoppingOffdiagonal
"""
struct EmbedHoppingOffdiagonal{C<:Number}
  amplitude ::C
  i ::Int64
  j ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
end


"""
    InteractionDiagonal
"""
struct EmbedInteractionDiagonal{R<:Real}
  amplitude ::R
  i ::Int64
  j ::Int64
  ri ::CarteCoord
  rj ::CarteCoord
end


"""
    InteractionOffdiagonal
"""
struct EmbedInteractionOffdiagonal{C<:Number}
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
    embed
"""
function embed(uc ::UnitCell{O},
               hopspec ::SpecHoppingDiagonal{R}) where {O, R<:Real}
  i = hopspec.i
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  return EmbedHoppingDiagonal{R}(hopspec.amplitude, i, ri)
end


"""
    embed
"""
function embed(uc ::UnitCell{O},
               hopspec ::SpecHoppingOffdiagonal{C}) where {O, C<:Number}
  (i, j) = (hopspec.i, hopspec.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)

  v = hopspec.amplitude
  @assert(i <= j, "ordering of i and j should have been taken care of by Spec")
  return EmbedHoppingOffdiagonal{C}(v, i, j, ri, rj)
end


"""
    embed
"""
function embed(uc::UnitCell{O},
               hopspec::SpecInteractionDiagonal{R}) where {O, R<:Real}
  (i, j) = (hopspec.i, hopspec.j)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)
  v = hopspec.amplitude
  @assert(i <= j, "ordering of i and j should have been taken care of by Spec")
  return EmbedInteractionDiagonal{R}(v, i, j, ri, rj)
end


"""
    embed
"""
function embed(uc::UnitCell{O},
               hopspec::SpecInteractionOffdiagonal{C}) where {O, C<:Number}
  (i, j) = (hopspec.i, hopspec.j)
  (k, l) = (hopspec.k, hopspec.l)
  ri = fract2carte(uc, getorbitalcoord(uc, hopspec.i) + hopspec.Ri)
  rj = fract2carte(uc, getorbitalcoord(uc, hopspec.j) + hopspec.Rj)
  rk = fract2carte(uc, getorbitalcoord(uc, hopspec.k) + hopspec.Rk)
  rl = fract2carte(uc, getorbitalcoord(uc, hopspec.l) + hopspec.Rl)
  v = hopspec.amplitude
  @assert(i <= j && k <= l && i <= k)
  return EmbedInteractionOffdiagonal{C}(v, i, j, k, l, ri, rj, rk, rl)
end


"""
    Hopping
"""
const EmbedHopping = Union{EmbedHoppingDiagonal, EmbedHoppingOffdiagonal}


"""
    Interaction
"""
const EmbedInteraction = Union{EmbedInteractionDiagonal, EmbedInteractionOffdiagonal}


"""
    Hamiltonian
"""
mutable struct EmbedHamiltonian{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{EmbedHopping}
  interactions ::Vector{EmbedInteraction}
end

"""
    Hamiltonian
"""
function EmbedHamiltonian(unitcell ::UnitCell{O},
                          hoppings ::AbstractVector{EmbedHopping},
                          interactions ::AbstractVector{EmbedInteraction}) where {O}
  return EmbedHamiltonian{O}(unitcell, hoppings, interactions)
end


"""
    Hamiltonian
"""
function embed(spec ::SpecHamiltonian{O}) where {O}
  unitcell = spec.unitcell
  hoppings = [embed(unitcell, hop) for hop in spec.hoppings]
  interactions = [embed(unitcell, inter) for inter in spec.interactions]
  EmbedHamiltonian{O}(unitcell, hoppings, interactions)
end


end # Embed
