module Realization

#=
type Term{S <: Number, N}
  amplitude ::S
  indices ::NTuple{N, Int64}
  positions ::NTuple{N, FractCoord}
end
=#

type HoppingDiagonal
  amplitude ::Real
  i ::Int64
  ri ::FractCoord
  function HoppingDiagonal(uc ::UnitCell, hopspec ::Spec.HoppingDiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    return new(hopspec.amplitude, i, ri + hopspec.Ri)
  end
end

type HoppingOffdiagonal{C <: Coord}
  amplitude ::Complex
  i, j ::Int64
  ri, rj ::FractCoord
  function HoppingOffdiagonal(uc ::UnitCell, hopspec ::Spec.HoppingOffdiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    (j, rj) = getorbitalindexcoord(uc, hopspec.j)
    return new(hopspec.amplitude, i, j, ri + hopspec.Ri, rj + hopspec.Rj)
  end
end

type InteractionDiagonal{C <: Coord}
  amplitude ::Complex
  i, j ::Int64
  ri, rj ::FractCoord
  function InteractionDiagonal(uc ::UnitCell, hopspec ::Spec.InteractionDiagonal)
    (i, ri) = getorbitalindexcoord(uc, hopspec.i)
    (j, rj) = getorbitalindexcoord(uc, hopspec.j)
    return new(hopspec.amplitude, i, j, ri + hopspec.Ri, rj + hopspec.Rj)
  end
end

type InteractionOffdiagonal{C <: Coord}
  amplitude ::Complex
  i, j, k, l ::Int64
  ri, rj, rk, rl ::FractCoord
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


end # Realization