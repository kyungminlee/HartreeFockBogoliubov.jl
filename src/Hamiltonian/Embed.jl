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
  if i > j
    (i,j) = (j,i)
    (ri, rj) = (rj, ri)
    v = conj(v)
  end
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
  if i > j
    (i,j) = (j,i)
    (ri, rj) = (rj, ri)
  end
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
  if i > j
    (i, j) = (j, i)
    (ri, rj) = (rj, ri)
    v = -v
  end

  if k > l
    (k, l) = (l, k)
    (rk, rl) = (rl, rk)
    v = -v
  end

  if i > k
    (i, j, k, l) = (k, l, i, j)
    (ri, rj, rk, rl) = (rk, rl, ri, rj)
    v = conj(v)
  end

  return InteractionOffdiagonal(v, i, j, k, l, ri, rj, rk, rl)
end


"""
    embed
"""
function embed(uc ::UnitCell, hopspec ::Spec.HoppingDiagonal)
  return HoppingDiagonal(uc, hopspec)
end


"""
    embed
"""
function embed(uc ::UnitCell, hopspec ::Spec.HoppingOffdiagonal)
  return HoppingOffdiagonal(uc, hopspec)
end


"""
    embed
"""
function embed(uc ::UnitCell, hopspec ::Spec.InteractionDiagonal)
  return InteractionDiagonal(uc, hopspec)
end


"""
    embed
"""
function embed(uc ::UnitCell, hopspec ::Spec.InteractionOffdiagonal)
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
type Hamiltonian
  unitcell ::UnitCell
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
end


"""
    Hamiltonian
"""
function Hamiltonian(spec ::Spec.Hamiltonian)
  unitcell = spec.unitcell
  hoppings = [embed(unitcell, hop) for hop in spec.hoppings]
  interactions = [embed(unitcell, inter) for inter in spec.interactions]
  Hamiltonian(unitcell, hoppings, interactions)
end

end # Embed
