module Spec

using ..Lattice


immutable HoppingDiagonal{T}
  amplitude ::Real
  i ::T
  Ri ::Vector{Int64}
end


immutable HoppingOffdiagonal{T}
  amplitude ::Number
  i ::T
  j ::T
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
  i ::T
  j ::T
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
  i ::T
  j ::T
  k ::T
  l ::T
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

"""
    hoppingbycarte
"""
function hoppingbycarte{T}(uc ::UnitCell{T},
                          amplitude ::Real,
                          i ::T,
                          ri ::CarteCoord;
                          tol=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  return HoppingDiagonal(amplitude, i, Ri)
end


function hoppingbycarte{T}(uc ::UnitCell{T},
                           amplitude ::Number,
                           i ::T, j ::T,
                           ri ::CarteCoord, rj ::CarteCoord;
                           tol=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return HoppingOffdiagonal{T}(amplitude, i, j, Ri, Rj)
end


function interactionbycarte{T}(uc ::UnitCell{T},
                              amplitude ::Real,
                              i ::T, j ::T,
                              ri ::CarteCoord, rj ::CarteCoord;
                              tol=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return InteractionDiagonal{T}(amplitude, i, j, Ri, Rj)
end


function interactionbycarte{T}(uc ::UnitCell{T},
                              amplitude ::Number,
                              i ::T, j ::T,
                              k ::T, l ::T,
                              ri ::CarteCoord, rj ::CarteCoord,
                              rk ::CarteCoord, rl ::CarteCoord;
                              tol=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  Rk = whichunitcell(uc, k, rk; tol=tol)
  Rl = whichunitcell(uc, l, rl; tol=tol)
  return InteractionOffdiagonal{T}(amplitude, i, j, k, l, Ri, Rj, Rk, Rl)
end


function addhopping!(hamiltonian ::Hamiltonian, hopping::HoppingDiagonal)
  @assert(haskey(hamiltonian.unitcell.orbitalindices, hopping.i),
          "orbital $(hopping.i) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


function addhopping!(hamiltonian ::Hamiltonian, hopping::HoppingOffdiagonal)
  @assert(haskey(hamiltonian.unitcell.orbitalindices, hopping.i),
          "orbital $(hopping.i) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, hopping.j),
          "orbital $(hopping.j) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


function addinteraction!(hamiltonian ::Hamiltonian, interaction::InteractionDiagonal)
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.i),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.j),
          "orbital $(interaction.j) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end


function addinteraction!(hamiltonian ::Hamiltonian, interaction::InteractionOffdiagonal)
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.i),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.j),
          "orbital $(interaction.j) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.k),
          "orbital $(interaction.k) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.l),
          "orbital $(interaction.l) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end


end # Spec
