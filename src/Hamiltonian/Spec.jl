"""
    Submodule `Spec`

"""
module Spec

using ..Lattice

export SpecHoppingDiagonal,
       SpecHoppingOffdiagonal,
       SpecInteractionDiagonal,
       SpecInteractionOffdiagonal
export SpecHopping, SpecInteraction
export SpecHamiltonian
export hoppingbycarte, interactionbycarte
export addhopping!, addinteraction!

"""
    SpecHoppingDiagonal{T}

Represents
```math
  t c_{i}^{*} c_{i}
```
  # Members
  * `amplitude ::Real`
  * `i ::Int64`: name of orbital
  * `Ri ::Vector{Int64}`: which unit cell? (indexed by a1, and a2)
"""
struct SpecHoppingDiagonal{R<:Real}
  amplitude ::R
  i ::Int64
  Ri ::Vector{Int64}
  function SpecHoppingDiagonal{R}(v::R,
                                  i::Integer,
                                  Ri::AbstractVector{Int64}) where {R<:Real}
    return new{R}(v, i, Ri)
  end
end

function SpecHoppingDiagonal(v::R,
                             i::Integer,
                             Ri::AbstractVector{Int64}) where {R<:Real}
  return SpecHoppingDiagonal{R}(v, i, Ri)
end


"""
    SpecHoppingOffdiagonal{T}

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
struct SpecHoppingOffdiagonal{C<:Number}
  amplitude ::C
  i ::Int64
  j ::Int64
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  function SpecHoppingOffdiagonal{C}(v::C,
                                     i::Integer,
                                     j::Integer,
                                     Ri::AbstractVector{Int64},
                                     Rj::AbstractVector{Int64}) where {C<:Number}
    @assert(length(Ri) == length(Rj))
    if i > j
      (i,j) = (j,i)
      (Ri, Rj) = (Rj, Ri)
      v = conj(v)
    end
    return new{C}(v, i, j, Ri, Rj)
  end
end


function SpecHoppingOffdiagonal(v::C,
                                i::Integer,
                                j::Integer,
                                Ri::AbstractVector{Int64},
                                Rj::AbstractVector{Int64}) where {C<:Number}
  return SpecHoppingOffdiagonal{C}(v, i, j, Ri, Rj)
end


"""
    SpecInteractionDiagonal{T}

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
struct SpecInteractionDiagonal{R<:Real}
  amplitude ::R
  i ::Int64
  j ::Int64
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  function SpecInteractionDiagonal{R}(v::R,
                                      i::Integer,
                                      j::Integer,
                                      Ri::AbstractVector{Int64},
                                      Rj::AbstractVector{Int64}) where {R<:Real}
    @assert(length(Ri) == length(Rj))
    @assert(!(i == j && Ri == Rj))
    if i > j
      (i,j) = (j,i)
      (Ri, Rj) = (Rj, Ri)
      v = conj(v)
    end
    return new{R}(v, i, j, Ri, Rj)
  end
end


function SpecInteractionDiagonal(v::R,
                                 i::Integer,
                                 j::Integer,
                                 Ri::AbstractVector{Int64},
                                 Rj::AbstractVector{Int64}) where {R<:Real}
  return SpecInteractionDiagonal{R}(v, i, j, Ri, Rj)
end


"""
    SpecInteractionOffdiagonal{T}

    i < j, k < l, i < k or (i == k and j < l)

Represents
```math
   U     c_{i}^{*} c_{j}^{*} c_{l} c_{k}
 + U^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}
```

Only keep the first term (and require i < j, k < l, i <= k)

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
struct SpecInteractionOffdiagonal{C<:Number}
  amplitude ::C
  i ::Int64
  j ::Int64
  k ::Int64
  l ::Int64
  Ri ::Vector{Int64}
  Rj ::Vector{Int64}
  Rk ::Vector{Int64}
  Rl ::Vector{Int64}
  function SpecInteractionOffdiagonal{C}(v::C,
                                         i::Integer,
                                         j::Integer,
                                         k::Integer,
                                         l::Integer,
                                         Ri::AbstractVector{Int64},
                                         Rj::AbstractVector{Int64},
                                         Rk::AbstractVector{Int64},
                                         Rl::AbstractVector{Int64}) where {C<:Number}
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
    return new{C}(v, i, j, k, l, Ri, Rj, Rk, Rl)
  end
end

function SpecInteractionOffdiagonal(v::C,
                                    i::Integer,
                                    j::Integer,
                                    k::Integer,
                                    l::Integer,
                                    Ri::AbstractVector{Int64},
                                    Rj::AbstractVector{Int64},
                                    Rk::AbstractVector{Int64},
                                    Rl::AbstractVector{Int64}) where {C<:Number}
  return SpecInteractionOffdiagonal{C}(v, i, j, k, l, Ri, Rj, Rk, Rl)
end


"""
    Hopping
"""
const SpecHopping = Union{SpecHoppingDiagonal, SpecHoppingOffdiagonal}


"""
    Interaction
"""
const SpecInteraction = Union{SpecInteractionDiagonal, SpecInteractionOffdiagonal}


"""
    Hamiltonian

  # Members
  * `unitcell ::UnitCell`
  * `hoppings ::Vector{Hopping}`
  * `interactions ::Vector{Interaction}`
"""
mutable struct SpecHamiltonian{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{SpecHopping}
  interactions ::Vector{SpecInteraction}
end


"""
    Hamiltonian

  Create an empty Hamiltonian

  # Arguments
  * `unitcell ::UnitCell`
"""
SpecHamiltonian(unitcell::UnitCell{O}) where {O} = SpecHamiltonian{O}(unitcell, [], [])


"""
    hoppingbycarte{T}

  # Arguments
  * `uc ::UnitCell{T}`
  * `amplitude ::Real`
  * `i ::T`
  * `ri ::CarteCoord`
  * `tol ::Real` : Optional. Defaults to `sqrt(eps(Float64))`
"""
function hoppingbycarte(uc ::UnitCell{O},
                        amplitude ::R,
                        i ::O,
                        ri ::CarteCoord;
                        tol::Real=sqrt(eps(Float64))) where {O, R<:Real}
  idx = getorbitalindex(uc, i)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  return SpecHoppingDiagonal{R}(amplitude, idx, Ri)
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
function hoppingbycarte(uc ::UnitCell{O},
                        amplitude ::C,
                        i ::O, j ::O,
                        ri ::CarteCoord, rj ::CarteCoord;
                        tol::Real=sqrt(eps(Float64))) where {O, C<:Number}
  idx = getorbitalindex(uc, i)
  jdx = getorbitalindex(uc, j)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return SpecHoppingOffdiagonal{C}(amplitude, idx, jdx, Ri, Rj)
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
function interactionbycarte(uc ::UnitCell{O},
                            amplitude ::R,
                            i ::O, j ::O,
                            ri ::CarteCoord, rj ::CarteCoord;
                            tol::Real=sqrt(eps(Float64))) where {O, R<:Real}
  idx = getorbitalindex(uc, i)
  jdx = getorbitalindex(uc, j)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  return SpecInteractionDiagonal{R}(amplitude, idx, jdx, Ri, Rj)
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
function interactionbycarte(uc ::UnitCell{O},
                            amplitude ::C,
                            i ::O, j ::O,
                            k ::O, l ::O,
                            ri ::CarteCoord, rj ::CarteCoord,
                            rk ::CarteCoord, rl ::CarteCoord;
                            tol=sqrt(eps(Float64))) where {O, C<:Number}
  idx = getorbitalindex(uc, i)
  jdx = getorbitalindex(uc, j)
  kdx = getorbitalindex(uc, k)
  ldx = getorbitalindex(uc, l)
  Ri = whichunitcell(uc, i, ri; tol=tol)
  Rj = whichunitcell(uc, j, rj; tol=tol)
  Rk = whichunitcell(uc, k, rk; tol=tol)
  Rl = whichunitcell(uc, l, rl; tol=tol)
  return SpecInteractionOffdiagonal{C}(amplitude, idx, jdx, kdx, ldx, Ri, Rj, Rk, Rl)
end


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingDiagonal`
"""
function addhopping!(hamiltonian ::SpecHamiltonian{O},
                     hopping ::SpecHoppingDiagonal{R};
                     tol=eps(Float64)) where {O, R<:Real}
  @assert(1 <= hopping.i <= numorbital(hamiltonian.unitcell),
          "orbital $(hopping.i) not defined in unit cell")
  if abs(hopping.amplitude) > tol
    push!(hamiltonian.hoppings, hopping)
  end
end


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingOffdiagonal`
"""
function addhopping!(hamiltonian ::SpecHamiltonian{O},
                     hopping ::SpecHoppingOffdiagonal{C};
                     tol=eps(Float64)) where {O, C<:Number}
  @assert(1 <= hopping.j <= numorbital(hamiltonian.unitcell),
          "orbital $(hopping.i) not defined in unit cell")
  @assert(1 <= hopping.j <= numorbital(hamiltonian.unitcell),
          "orbital $(hopping.j) not defined in unit cell")
  if abs(hopping.amplitude) > tol
    push!(hamiltonian.hoppings, hopping)
  end
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionDiagonal`
"""
function addinteraction!(hamiltonian ::SpecHamiltonian{O},
                         interaction ::SpecInteractionDiagonal{R};
                         tol=eps(Float64)) where {O, R<:Real}
  @assert(1 <= interaction.i <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(1 <= interaction.j <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.j) not defined in unit cell")
  if abs(interaction.amplitude) > tol
    push!(hamiltonian.interactions, interaction)
  end
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionOffdiagonal`
"""
function addinteraction!(hamiltonian ::SpecHamiltonian{O},
                         interaction::SpecInteractionOffdiagonal{C};
                         tol=eps(Float64)) where {O, C<:Number}
  @assert(1 <= interaction.i <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.i) not defined in unit cell")
  @assert(1 <= interaction.j <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.j) not defined in unit cell")
  @assert(1 <= interaction.k <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.k) not defined in unit cell")
  @assert(1 <= interaction.l <= numorbital(hamiltonian.unitcell),
          "orbital $(interaction.l) not defined in unit cell")
  if abs(interaction.amplitude) > tol
    push!(hamiltonian.interactions, interaction)
  end
end

end # Spec
