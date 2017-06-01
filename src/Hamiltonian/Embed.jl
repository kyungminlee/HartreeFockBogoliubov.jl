"""
Submodule `Embed`

From `Spec` for Hamiltonian in periodic system, create an `Embed`


"""
module Embed

#import HartreeFockBogoliubov.PeriodicSystem: CarteCoord, FractCoord, UnitCell
using ..Lattice
import ..Spec

#import HartreeFockBogoliubov.PeriodicSystem
#=
type Term{S <: Number, N}
  amplitude ::S
  indices ::NTuple{N, Int64}
  positions ::NTuple{N, FractCoord}
end
=#

"""
    HoppingDiagonal
"""
immutable HoppingDiagonal
  amplitude ::Real
  i ::Int64
  ri ::FractCoord
end

function HoppingDiagonal(uc ::UnitCell, hopspec ::Spec.HoppingDiagonal)
  (i, ri) = getorbitalindexcoord(uc, hopspec.i)
  return HoppingDiagonal(hopspec.amplitude, i, ri + hopspec.Ri)
end


"""
    HoppingOffdiagonal
"""
immutable HoppingOffdiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  ri ::FractCoord
  rj ::FractCoord
end

function HoppingOffdiagonal(uc ::UnitCell, hopspec ::Spec.HoppingOffdiagonal)
  (i, ri) = getorbitalindexcoord(uc, hopspec.i)
  (j, rj) = getorbitalindexcoord(uc, hopspec.j)
  (Ri, Rj) = (hopspec.Ri, hopspec.Rj)
  v = hopspec.amplitude
  if i > j
    (i,j) = (j,i)
    (ri, rj) = (rj, ri)
    (Ri, Rj) = (Rj, Ri)
    v = conj(v)
  end
  return HoppingOffdiagonal(v, i, j, ri + Ri, rj + Rj)
end

immutable InteractionDiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  ri ::FractCoord
  rj ::FractCoord
end

function InteractionDiagonal(uc ::UnitCell, hopspec ::Spec.InteractionDiagonal)
  (i, ri) = getorbitalindexcoord(uc, hopspec.i)
  (j, rj) = getorbitalindexcoord(uc, hopspec.j)
  (Ri, Rj) = (hopspec.Ri, hopspec.Rj)
  v = hopspec.amplitude
  if i > j
    (i,j) = (j,i)
    (ri, rj) = (rj, ri)
    (Ri, Rj) = (Rj, Ri)
  end
  return InteractionDiagonal(v, i, j, ri + Ri, rj + Rj)
end


immutable InteractionOffdiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  k ::Int64
  l ::Int64
  ri ::FractCoord
  rj ::FractCoord
  rk ::FractCoord
  rl ::FractCoord
  function InteractionOffdiagonal(uc ::UnitCell, hopspec ::Spec.InteractionOffdiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    (j, rj) = getorbitalindexcoord(uc, hopspec.j)
    (k, rk) = getorbitalindexcoord(uc, hopspec.k)
    (l, rl) = getorbitalindexcoord(uc, hopspec.l)
    (Ri, Rj) = (hopspec.Ri, hopspec.Rj)
    (Rk, Rl) = (hopspec.Rk, hopspec.Rl)
    v = hopspec.amplitude
    if i > j
      (i, j) = (j, i)
      (ri, rj) = (rj, ri)
      (Ri, Rj) = (Rj, Ri)
      v = -v
    end

    if k > l
      (k, l) = (l, k)
      (rk, rl) = (rl, rk)
      (Rk, Rl) = (Rl, Rk)
      v = -v
    end

    if i > k
      (i, j, k, l) = (k, l, i, j)
      (ri, rj, rk, rl) = (rk, rl, ri, rj)
      (Ri, Rj, Rk, Rl) = (Rk, Rl, Ri, Rj)
      v = conj(v)
    end

    return InteractionOffdiagonal(v,
               i, j, k, l,
               ri + Ri, rj + Rj, rk + Rk, rl + Rl)
  end
end



function embed(uc ::UnitCell, hopspec ::Spec.HoppingDiagonal)
  return HoppingDiagonal(uc, hopspec)
end

function embed(uc ::UnitCell, hopspec ::Spec.HoppingOffdiagonal)
  return HoppingOffdiagonal(uc, hopspec)
end

function embed(uc ::UnitCell, hopspec ::Spec.InteractionDiagonal)
  return InteractionDiagonal(uc, hopspec)
end

function embed(uc ::UnitCell, hopspec ::Spec.InteractionOffdiagonal)
  return InteractionOffdiagonal(uc, hopspec)
end


typealias Hopping Union{HoppingDiagonal, HoppingOffdiagonal}

typealias Interaction Union{InteractionDiagonal, InteractionOffdiagonal}


type Hamiltonian
  unitcell ::UnitCell
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
end


end # Embed
