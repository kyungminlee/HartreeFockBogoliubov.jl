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
type HoppingDiagonal
  amplitude ::Real
  i ::Int64
  ri ::FractCoord
  function HoppingDiagonal(uc ::UnitCell, hopspec ::Spec.HoppingDiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    return new(hopspec.amplitude, i, ri + hopspec.Ri)
  end
end


"""
    HoppingOffdiagonal
"""
type HoppingOffdiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  ri ::FractCoord
  rj ::FractCoord
  function HoppingOffdiagonal(uc ::UnitCell, hopspec ::Spec.HoppingOffdiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    (j, rj) = getorbitalindexcoord(uc, hopspec.j)
    return new(hopspec.amplitude, i, j, ri + hopspec.Ri, rj + hopspec.Rj)
  end
end



type InteractionDiagonal
  amplitude ::Complex
  i ::Int64
  j ::Int64
  ri ::FractCoord
  rj ::FractCoord
  function InteractionDiagonal(uc ::UnitCell, hopspec ::Spec.InteractionDiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    (j, rj) = getorbitalindexcoord(uc, hopspec.j)
    return new(hopspec.amplitude, i, j, ri + hopspec.Ri, rj + hopspec.Rj)
  end
end



type InteractionOffdiagonal
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
    return new(hopspec.amplitude,
               i, j, k, l,
               ri + hopspec.Ri,
               rj + hopspec.Rj,
               rk + hopspec.Rk,
               rl + hopspec.Rl,
               )
  end
end


typealias Hopping Union{HoppingDiagonal, HoppingOffdiagonal}


typealias Interaction Union{InteractionDiagonal, InteractionOffdiagonal}


type Hamiltonian
  unitcell ::UnitCell
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
end


end # Embed