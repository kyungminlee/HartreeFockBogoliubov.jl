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
immutable HoppingDiagonal
  amplitude ::Real
  i ::Int64
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
immutable HoppingOffdiagonal
  amplitude ::Number
  i ::Int64
  j ::Int64
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  function HoppingOffdiagonal(v::Number, i::Int64, j::Int64, Ri::Vector{Int64}, Rj::Vector{Int64})
    @assert(length(Ri) == length(Rj))
    if i > j
      (i,j) = (j,i)
      (Ri, Rj) = (Rj, Ri)
      v = conj(v)
    end
    return new(v, i, j, Ri, Rj)
  end
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
immutable InteractionDiagonal
  amplitude ::Real
  i ::Int64
  j ::Int64
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  function InteractionDiagonal(v::Real,
                               i::Int64,
                               j::Int64,
                               Ri::AbstractVector{Int64},
                               Rj::AbstractVector{Int64})
    @assert(length(Ri) == length(Rj))
    @assert(!(i == j && Ri == Rj))
    if i > j
      (i,j) = (j,i)
      (Ri, Rj) = (Rj, Ri)
      v = conj(v)
    end
    return new(v, i, j, Ri, Rj)
  end
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
immutable InteractionOffdiagonal
  amplitude ::Number
  i ::Int64
  j ::Int64
  k ::Int64
  l ::Int64
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  Rk ::Vector{Int64}
  Rl ::Vector{Int64}
  function InteractionOffdiagonal(v::Number,
                                  i::Int64,
                                  j::Int64,
                                  k::Int64,
                                  l::Int64,
                                  Ri::AbstractVector{Int64},
                                  Rj::AbstractVector{Int64},
                                  Rk::AbstractVector{Int64},
                                  Rl::AbstractVector{Int64})
    @assert(length(Ri) == length(Rj) == length(Rk) == length(Rl))
    @assert(i != j)
    @assert(k != l)
    if i > j
      (i, j) = (j, i)
      (Ri, Rj) = (Rj, Ri)
      v = -v
    end

    if k > l
      (k, l) = (l, k)
      (Rk, Rl) = (Rl, Rk)
      v = -v
    end

    if i > k
      (i, j, k, l) = (k, l, i, j)
      (Ri, Rj, Rk, Rl) = (Rk, Rl, Ri, Rj)
      v = conj(v)
    end
    return new(v, i, j, k, l, Ri, Rj, Rk, Rl)
  end
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

  # Members
  * `unitcell ::UnitCell`
  * `hoppings ::Vector{Hopping}`
  * `interactions ::Vector{Interaction}`
"""
type Hamiltonian{T}
  unitcell ::UnitCell{T}
  hoppings ::Vector{Hopping}
  interactions ::Vector{Interaction}
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
  idx = getorbitalindex(uc, i)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  return HoppingDiagonal(amplitude, idx, Ri)
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
  idx = getorbitalindex(uc, i)
  jdx = getorbitalindex(uc, j)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return HoppingOffdiagonal(amplitude, idx, jdx, Ri, Rj)
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
  idx = getorbitalindex(uc, i)
  jdx = getorbitalindex(uc, j)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return InteractionDiagonal(amplitude, idx, jdx, Ri, Rj)
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
  idx = getorbitalindex(uc, i)
  jdx = getorbitalindex(uc, j)
  kdx = getorbitalindex(uc, k)
  ldx = getorbitalindex(uc, l)
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
function addhopping!{T}(hamiltonian ::Hamiltonian{T}, hopping ::HoppingDiagonal; tol=eps(Float64))
  @assert(1 <= hopping.i <= numorbital(hamiltonian.unitcell),
          "orbital $(hopping.i) not defined in unit cell")
  #if abs(hopping.amplitude) > tol
    push!(hamiltonian.hoppings, hopping)
  #end
end


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingOffdiagonal`
"""
function addhopping!{T}(hamiltonian ::Hamiltonian{T}, hopping ::HoppingOffdiagonal; tol=eps(Float64))
  @assert(1 <= hopping.j <= numorbital(hamiltonian.unitcell),
          "orbital $(hopping.i) not defined in unit cell")
  @assert(1 <= hopping.j <= numorbital(hamiltonian.unitcell),
          "orbital $(hopping.j) not defined in unit cell")
  #if abs(hopping.amplitude) > tol
    push!(hamiltonian.hoppings, hopping)
  #end
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionDiagonal`
"""
function addinteraction!{T}(hamiltonian ::Hamiltonian{T}, interaction ::InteractionDiagonal; tol=eps(Float64))
  @assert(1 <= interaction.i <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(1 <= interaction.j <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.j) not defined in unit cell")
  #if abs(interaction.amplitude) > tol
    push!(hamiltonian.interactions, interaction)
  #end
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionOffdiagonal`
"""
function addinteraction!{T}(hamiltonian ::Hamiltonian{T}, interaction::InteractionOffdiagonal; tol=eps(Float64))
  @assert(1 <= interaction.i <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(1 <= interaction.j <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.j) not defined in unit cell")
  @assert(1 <= interaction.k <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.k) not defined in unit cell")
  @assert(1 <= interaction.l <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.l) not defined in unit cell")
  #if abs(interaction.amplitude) > tol
    push!(hamiltonian.interactions, interaction)
  #end
end


end # Spec
