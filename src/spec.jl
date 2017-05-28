module Spec

#type Term{S <: Number, N}
#  amplitude ::S
#  orbital ::NTuple{String, N}
#  unitcelldisplacement ::Vector{Int64}
#end

type HoppingDiagonal
  amplitude ::Real
  i ::String
  Ri ::Vector{Int64}
end

type HoppingOffdiagonal
  amplitude ::Number
  i, j ::String
  Ri, Rj ::Vector{Int64}
end

type InteractionDiagonal
  amplitude ::Real
  i, j ::String
  Ri, Rj ::Vector{Int64}
end

type InteractionOffdiagonal
  amplitude ::Number
  i, j, k, l ::String
  Ri, Rj, Rk, Rl ::Vector{Int64}
end

typealias Hopping Union{HoppingDiagonal, HoppingOffdiagonal}
typealias Interaction Union{InteractionDiagonal, InteractionOffdiagonal}

type Hamiltonian
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
end

end Spec
