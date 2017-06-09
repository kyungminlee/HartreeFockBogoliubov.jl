"""
    Submodule `Spec`

"""
module Spec

using ..Lattice

"""
    HoppingDiagonal{T}

Represents
```math
  t c_{i}^{*} c_{i}
```
  # Members
  * `amplitude ::Real`
  * `i ::T`
  * `Ri ::Vector{Int64}`
"""
immutable HoppingDiagonal{T}
  amplitude ::Real
  i ::T
  Ri ::Vector{Int64}
end


"""
    HoppingOffdiagonal{T}

Represents
```math
  t c_{i}^{*} c_{j} + t^* c_{j}^{*} c_{i}
```

  # Members
  * `amplitude ::Number`
  * `i ::T`
  * `j ::T`
  * `Ri ::Vector{Int64}`
  * `Rj ::Vector{Int64}`
"""
immutable HoppingOffdiagonal{T}
  amplitude ::Number
  i ::T
  j ::T
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
end


"""
    InteractionDiagonal{T}

Represents
```math
    U c_{i}^{*} c_{j}^{*} c_{j} c_{i}
```

  # Members
  * `amplitude ::Real`
  * `i ::T`
  * `j ::T`
  * `Ri ::Vector{Int64}`
  * `Rj ::Vector{Int64}`
"""
immutable InteractionDiagonal{T}
  amplitude ::Real
  i ::T
  j ::T
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
end


"""
    InteractionOffdiagonal{T}

    i < j, k < l, i < k or (i == k and j < l)

Represents
```math
   U     c_{i}^{*} c_{j}^{*} c_{l} c_{k}
 + U^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}
```

  # Members
  * `amplitude ::Number`
  * `i ::T`
  * `j ::T`
  * `k ::T`
  * `l ::T`
  * `Ri ::Vector{Int64}`
  * `Rj ::Vector{Int64}`
  * `Rk ::Vector{Int64}`
  * `Rl ::Vector{Int64}`
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


"""
    Hopping
"""
const Hopping{T} = Union{HoppingDiagonal{T}, HoppingOffdiagonal{T}}


"""
    Interaction
"""
const Interaction{T} = Union{InteractionDiagonal{T}, InteractionOffdiagonal{T}}


"""
    Hamiltonian

  # Members
  * `unitcell ::UnitCell`
  * `hoppings ::Vector{Hopping}`
  * `interactions ::Vector{Interaction}`
"""
type Hamiltonian{T}
  unitcell ::UnitCell{T}
  hoppings ::Vector{Hopping{T}}
  interactions ::Vector{Interaction{T}}
end

"""
    Hamiltonian

  Create an empty Hamiltonian

  # Arguments
  * `unitcell ::UnitCell`
"""
Hamiltonian{T}(unitcell ::UnitCell{T}) = Hamiltonian{T}(unitcell, [], [])


"""
    hoppingbycarte{T}

  # Arguments
  * `uc ::UnitCell{T}`
  * `amplitude ::Real`
  * `i ::T`
  * `ri ::CarteCoord`
  * `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`
"""
function hoppingbycarte{T}(uc ::UnitCell{T},
                           amplitude ::Real,
                           i ::T,
                           ri ::CarteCoord;
                           tol::Real=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  return HoppingDiagonal(amplitude, i, Ri)
end


"""
    hoppingbycarte{T}

  # Arguments
  * `uc ::UnitCell{T}`
  * `amplitude ::Number`
  * `i ::T`
  * `j ::T`
  * `ri ::CarteCoord`
  * `rj ::CarteCoord`
  * `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`
"""
function hoppingbycarte{T}(uc ::UnitCell{T},
                           amplitude ::Number,
                           i ::T, j ::T,
                           ri ::CarteCoord, rj ::CarteCoord;
                           tol::Real=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return HoppingOffdiagonal{T}(amplitude, i, j, Ri, Rj)
end


"""
    interactionbycarte{T}

  # Arguments
    * `uc ::UnitCell{T}`
    * `amplitude ::Number`
    * `i ::T`
    * `j ::T`
    * `ri ::CarteCoord`
    * `rj ::CarteCoord`
    * `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`
"""
function interactionbycarte{T}(uc ::UnitCell{T},
                              amplitude ::Real,
                              i ::T, j ::T,
                              ri ::CarteCoord, rj ::CarteCoord;
                              tol::Real=sqrt(eps(Float64)))
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return InteractionDiagonal{T}(amplitude, i, j, Ri, Rj)
end


"""
    interactionbycarte{T}

  # Arguments
    * `uc ::UnitCell{T}`
    * `amplitude ::Number`
    * `i ::T`
    * `j ::T`
    * `k ::T`
    * `l ::T`
    * `ri ::CarteCoord`
    * `rj ::CarteCoord`
    * `rk ::CarteCoord`
    * `rl ::CarteCoord`
    * `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`
"""
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


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingDiagonal`
"""
function addhopping!{T}(hamiltonian ::Hamiltonian{T}, hopping ::HoppingDiagonal{T})
  @assert(haskey(hamiltonian.unitcell.orbitalindices, hopping.i),
          "orbital $(hopping.i) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingOffdiagonal`
"""
function addhopping!{T}(hamiltonian ::Hamiltonian{T}, hopping ::HoppingOffdiagonal{T})
  @assert(haskey(hamiltonian.unitcell.orbitalindices, hopping.i),
          "orbital $(hopping.i) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, hopping.j),
          "orbital $(hopping.j) not defined in unit cell")
  push!(hamiltonian.hoppings, hopping)
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionDiagonal`
"""
function addinteraction!{T}(hamiltonian ::Hamiltonian{T}, interaction ::InteractionDiagonal{T})
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.i),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(haskey(hamiltonian.unitcell.orbitalindices, interaction.j),
          "orbital $(interaction.j) not defined in unit cell")
  push!(hamiltonian.interactions, interaction)
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionOffdiagonal`
"""
function addinteraction!{T}(hamiltonian ::Hamiltonian{T}, interaction::InteractionOffdiagonal{T})
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
