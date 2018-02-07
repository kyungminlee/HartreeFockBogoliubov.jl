"""
    Submodule `Spec`

"""
module Spec

using ..Lattice

export HoppingDiagonal,
       HoppingOffdiagonal,
       InteractionDiagonal,
       InteractionOffdiagonal
export Hopping, Interaction
export FullHamiltonian
export hoppingbycarte, interactionbycarte
export addhopping!, addinteraction!
export islocal, localize

"""
    HoppingDiagonal{R<:Real}

Represents
```math
  t c_{i}^{*} c_{i}
```
  # Members
  * `amplitude ::R`
  * `i ::Int64`: name of orbital
  * `Ri ::Vector{Int64}`: which unit cell? (indexed by a1, and a2)
"""
struct HoppingDiagonal{R<:Real}
    amplitude ::R
    i ::Int64
    Ri ::Vector{Int64}
    function HoppingDiagonal{R}(v::R,
                                i::Integer,
                                Ri::AbstractVector{<:Integer}) where {R<:Real}
        return new{R}(v, i, Ri)
    end
end

"""

"""
function HoppingDiagonal(v::R,
                         i::Integer,
                         Ri::AbstractVector{<:Integer}) where {R<:Real}
    return HoppingDiagonal{R}(v, i, Ri)
end


"""
    HoppingOffdiagonal{C<:Number}

Represents
```math
  t c_{i}^{*} c_{j} + t^* c_{j}^{*} c_{i}
```

  # Members
  * `amplitude :: C`
  * `i ::T`
  * `j ::T`
  * `Ri ::Vector{Int64}`
  * `Rj ::Vector{Int64}`
"""
struct HoppingOffdiagonal{C<:Number}
    amplitude ::C
    i ::Int64
    j ::Int64
    Ri ::Vector{Int64}
    Rj ::Vector{Int64}
    function HoppingOffdiagonal{C}(v::C,
                                   i::Integer,
                                   j::Integer,
                                   Ri::AbstractVector{<:Integer},
                                   Rj::AbstractVector{<:Integer}) where {C<:Number}
        @assert(length(Ri) == length(Rj))
        if i > j
            (i,j) = (j,i)
            (Ri, Rj) = (Rj, Ri)
            v = conj(v)
        end
        return new{C}(v, i, j, Ri, Rj)
    end
end


function HoppingOffdiagonal(v::C,
                            i::Integer,
                            j::Integer,
                            Ri::AbstractVector{<:Integer},
                            Rj::AbstractVector{<:Integer}) where {C<:Number}
    return HoppingOffdiagonal{C}(v, i, j, Ri, Rj)
end


"""
    InteractionDiagonal{R<:Real}

Represents
```math
    U c_{i}^{*} c_{j}^{*} c_{j} c_{i}
```

  # Members
  * `amplitude ::R`
  * `i ::T`
  * `j ::T`
  * `Ri ::Vector{Int64}`
  * `Rj ::Vector{Int64}`
"""
struct InteractionDiagonal{R<:Real}
    amplitude ::R
    i ::Int64
    j ::Int64
    Ri ::Vector{Int64}
    Rj ::Vector{Int64}
    function InteractionDiagonal{R}(v::R,
                                    i::Integer,
                                    j::Integer,
                                    Ri::AbstractVector{<:Integer},
                                    Rj::AbstractVector{<:Integer}) where {R<:Real}
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


function InteractionDiagonal(v::R,
                             i::Integer,
                             j::Integer,
                             Ri::AbstractVector{<:Integer},
                             Rj::AbstractVector{<:Integer}) where {R<:Real}
    return InteractionDiagonal{R}(v, i, j, Ri, Rj)
end


"""
    InteractionOffdiagonal{C<:Number}

    i < j, k < l, i < k or (i == k and j < l)

Represents
```math
   U     c_{i}^{*} c_{j}^{*} c_{l} c_{k}
 + U^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}
```

Only keep the first term (and require i < j, k < l, i <= k)

  # Members
  * `amplitude ::C`
  * `i ::T`
  * `j ::T`
  * `k ::T`
  * `l ::T`
  * `Ri ::Vector{Int64}`
  * `Rj ::Vector{Int64}`
  * `Rk ::Vector{Int64}`
  * `Rl ::Vector{Int64}`
"""
struct InteractionOffdiagonal{C<:Number}
    amplitude ::C
    i ::Int64
    j ::Int64
    k ::Int64
    l ::Int64
    Ri ::Vector{Int64}
    Rj ::Vector{Int64}
    Rk ::Vector{Int64}
    Rl ::Vector{Int64}
    function InteractionOffdiagonal{C}(v::C,
                                       i::Integer,
                                       j::Integer,
                                       k::Integer,
                                       l::Integer,
                                       Ri::AbstractVector{<:Integer},
                                       Rj::AbstractVector{<:Integer},
                                       Rk::AbstractVector{<:Integer},
                                       Rl::AbstractVector{<:Integer}) where {C<:Number}
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

function InteractionOffdiagonal(v::C,
                                i::Integer,
                                j::Integer,
                                k::Integer,
                                l::Integer,
                                Ri::AbstractVector{<:Integer},
                                Rj::AbstractVector{<:Integer},
                                Rk::AbstractVector{<:Integer},
                                Rl::AbstractVector{<:Integer}) where {C<:Number}
    return InteractionOffdiagonal{C}(v, i, j, k, l, Ri, Rj, Rk, Rl)
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
    FullHamiltonian

  # Members
  * `unitcell ::UnitCell`
  * `hoppings ::Vector{Hopping}`
  * `interactions ::Vector{Interaction}`
"""
mutable struct FullHamiltonian{O}
    unitcell ::UnitCell{O}
    hoppings ::Vector{Hopping}
    interactions ::Vector{Interaction}
end


"""
    Hamiltonian

  Create an empty Hamiltonian

  # Arguments
  * `unitcell ::UnitCell`
"""
FullHamiltonian(unitcell::UnitCell{O}) where {O} = FullHamiltonian{O}(unitcell, [], [])


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
    return HoppingDiagonal{R}(amplitude, idx, Ri)
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
    return HoppingOffdiagonal{C}(amplitude, idx, jdx, Ri, Rj)
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
    return InteractionDiagonal{R}(amplitude, idx, jdx, Ri, Rj)
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
                            tol::Real=sqrt(eps(Float64))) where {O, C<:Number}
    idx = getorbitalindex(uc, i)
    jdx = getorbitalindex(uc, j)
    kdx = getorbitalindex(uc, k)
    ldx = getorbitalindex(uc, l)
    Ri = whichunitcell(uc, i, ri; tol=tol)
    Rj = whichunitcell(uc, j, rj; tol=tol)
    Rk = whichunitcell(uc, k, rk; tol=tol)
    Rl = whichunitcell(uc, l, rl; tol=tol)
    return InteractionOffdiagonal{C}(amplitude, idx, jdx, kdx, ldx, Ri, Rj, Rk, Rl)
end



function islocal(hopping::HoppingDiagonal{R}) where {R <: Real}
    return iszero(hopping.Ri)
end

function islocal(hopping::HoppingOffdiagonal{C}) where {C <: Number}
    return iszero(hopping.Ri)
end

function islocal(interaction::InteractionDiagonal{R}) where {R <: Real}
    return iszero(hopping.Ri)
end

function islocal(interaction::InteractionOffdiagonal{C}) where {C <: Number}
    return iszero(hopping.Ri)
end


function localize(hopping::HoppingDiagonal{R}) where {R <: Real}
    return HoppingDiagonal{R}(hopping.amplitude, hopping.i, zeros(hopping.Ri))
end

function localize(hopping::HoppingOffdiagonal{C}) where {C <: Number}
    return HoppingOffdiagonal{C}(hopping.amplitude,
                                 hopping.i,
                                 hopping.j,
                                 zeros(hopping.Ri),
                                 hopping.Rj - hopping.Ri)
end

function localize(interaction::InteractionDiagonal{R}) where {R <: Real}
    return InteractionDiagonal{R}(interaction.amplitude,
                                  interaction.i,
                                  interaction.j,
                                  zeros(interaction.Ri),
                                  interaction.Rj - interaction.Ri)
end

function localize(interaction::InteractionOffdiagonal{C}) where {C <: Number}
    return InteractionOffdiagonal{C}(interaction.amplitude,
                                     interaction.i,
                                     interaction.j,
                                     interaction.k,
                                     interaction.l,
                                     zeros(interaction.Ri),
                                     interaction.Rj - interaction.Ri,
                                     interaction.Rk - interaction.Ri,
                                     interaction.Rl - interaction.Ri)
end


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingDiagonal`
"""
function addhopping!(hamiltonian ::FullHamiltonian{O},
                     hopping ::HoppingDiagonal{R};
                     tol::Real=eps(Float64)) where {O, R<:Real}
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
function addhopping!(hamiltonian ::FullHamiltonian{O},
                     hopping ::HoppingOffdiagonal{C};
                     tol::Real=eps(Float64)) where {O, C<:Number}
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
function addinteraction!(hamiltonian ::FullHamiltonian{O},
                         interaction ::InteractionDiagonal{R};
                         tol::Real=eps(Float64)) where {O, R<:Real}
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
function addinteraction!(hamiltonian ::FullHamiltonian{O},
                         interaction::InteractionOffdiagonal{C};
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


import Base:isless

isless(hop1 ::HoppingDiagonal, hop2::HoppingOffdiagonal) = true
isless(hop1 ::HoppingOffdiagonal, hop2::HoppingDiagonal) = false

function simplify(hamiltonian::FullHamiltonian{O}) where {O}

    hopdiadict = Dict{Int64, Float64}()
    hopoffdict = Dict{Tuple{Int64, Int64, Vector{Int64}}, Complex128}()

    for hop in hoppings
        if isa(newhop, HoppingDiagonal)
            k = hop.i
            if haskey(hopdiadict, k)
                hopdiadict[k].amplitude += hop.amplitude
            else
                hopdiadict[k] = hop
            end
        else
            k = (hop.i, hop.j, hop.Rj - hop.Ri)
            if haskey(hopoffdict, k)
                hopoffdict[k].amplitude += hop.amplitude
            else
                hopoffdict[k] = hop
            end
        end
    end

    # TODO
    error("unimplemented")


    #newhamiltonian = FullHamiltonian{O}(hamiltonian.unitcell)
    #for hop in hamiltonian.hoppings
    #    addhopping!
    #end


end

end # Spec
