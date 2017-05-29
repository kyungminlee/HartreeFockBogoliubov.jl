module Spec

using ..Lattice

type Term{S <: Number, N}
  amplitude ::S
  subscript ::NTuple{N, AbstractString}
  latticedisplacement ::Vector{Int64}
end


type HoppingDiagonal
  amplitude ::Real
  i ::String
  Ri ::Vector{Int64}
end


type HoppingOffdiagonal
  amplitude ::Number
  i ::String
  j ::String
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
end

"""
i < j

Represents
```math
  U c_{i}^{*} c_{j}^{*} c_{j} c_{i}
```
"""
type InteractionDiagonal
  amplitude ::Real
  i ::String
  j ::String
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
end


"""
i < j, k < l, i < k or (i == k and j < l)

Represents 
```math
   U     c_{i}^{*} c_{j}^{*} c_{l} c_{k}
 + U^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}
```
"""
type InteractionOffdiagonal
  amplitude ::Number
  i ::String
  j ::String
  k ::String
  l ::String
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  Rk ::Vector{Int64}
  Rl ::Vector{Int64}
end


typealias Hopping Union{HoppingDiagonal, HoppingOffdiagonal}


typealias Interaction Union{InteractionDiagonal, InteractionOffdiagonal}


type Hamiltonian
  unitcell ::UnitCell
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
  Hamiltonian(unitcell ::UnitCell) = new(unitcell, [], [])
end


function addhopping!(hamiltonian ::Hamiltonian, hopping::HoppingDiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, hopping.i), "orbital $(hopping.i) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


function addhopping!(hamiltonian ::Hamiltonian, hopping::HoppingOffdiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, hopping.i), "orbital $(hopping.i) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, hopping.j), "orbital $(hopping.j) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


function addinteraction!(hamiltonian ::Hamiltonian, interaction::InteractionDiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.i), "orbital $(interaction.i) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.j), "orbital $(interaction.j) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end


function addinteraction!(hamiltonian ::Hamiltonian, interaction::InteractionOffdiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.i), "orbital $(interaction.i) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.j), "orbital $(interaction.j) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.k), "orbital $(interaction.k) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.l), "orbital $(interaction.l) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end



end # Spec
