module Spec

using ..Lattice

type Term{S <: Number, N}
  amplitude ::S
  subscript ::NTuple{N, AbstractString}
  latticedisplacement ::Vector{Int64}
end


immutable HoppingDiagonal{T}
  amplitude ::Real
  i ::String
  Ri ::Vector{Int64}
end


immutable HoppingOffdiagonal{T}
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
immutable InteractionDiagonal{T}
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
immutable InteractionOffdiagonal{T}
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


function hoppingbycarte(uc ::UnitCell, amplitude ::Real, i ::AbstractString, ri ::CarteCoord)
  Ri = whichunitcell(uc, i, ri)
  return HoppingDiagonal(amplitude, i, Ri)
end

function hoppingbycarte(uc ::UnitCell, amplitude ::Number,
                        i ::AbstractString, j ::AbstractString, 
                        ri ::CarteCoord, rj ::CarteCoord)
  Ri = whichunitcell(uc, i, ri)
  Rj = whichunitcell(uc, j, rj)
  return HoppingOffdiagonal(amplitude, i, j, Ri, Rj)
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
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, hopping.i),
          "orbital $(hopping.i) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


function addhopping!(hamiltonian ::Hamiltonian, hopping::HoppingOffdiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, hopping.i),
          "orbital $(hopping.i) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, hopping.j),
          "orbital $(hopping.j) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


function addinteraction!(hamiltonian ::Hamiltonian, interaction::InteractionDiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.i),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.j),
          "orbital $(interaction.j) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end


function addinteraction!(hamiltonian ::Hamiltonian, interaction::InteractionOffdiagonal)
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.i),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.j),
          "orbital $(interaction.j) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.k),
          "orbital $(interaction.k) not defined in unit cell")
  @assert(!haskey(hamiltonian.unitcell.orbitalindices, interaction.l),
          "orbital $(interaction.l) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end



end # Spec
