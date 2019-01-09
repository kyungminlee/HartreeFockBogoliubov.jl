"""
    HoppingDiagonal{R<:Real}

Represents
```math
  t c_{i}^{*} c_{i}
```

  # Members
  * `amplitude ::R`
  * `i ::Int`: index of orbital
  * `Ri ::Vector{Int}`: which unit cell? (indexed by a1, and a2)
"""
struct HoppingDiagonal{R<:Real}
    amplitude ::R
    i ::Int
    Ri ::Vector{Int}
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

`t`, `i`, `j` and the unitcell-coordinates `Ri` and `Rj` are stored.
Require that `(i, Ri) <= (j, Rj)`


# Members
* `amplitude ::C`
* `i, j ::T`
* `Ri, Rj ::Vector{Int}`
"""
struct HoppingOffdiagonal{C<:Number}
    amplitude ::C
    i ::Int
    j ::Int
    Ri ::Vector{Int}
    Rj ::Vector{Int}
    function HoppingOffdiagonal{C}(v::C,
                                   i::Integer,
                                   j::Integer,
                                   Ri::AbstractVector{<:Integer},
                                   Rj::AbstractVector{<:Integer}) where {C<:Number}
        if length(Ri) != length(Rj)
            throw(ArgumentError("Ri and Rj must be of same length"))
        elseif (i, Ri) == (j, Rj)
            throw(ArgumentError("not an off-diagonal element"))
        end

        if (i, Ri) > (j, Rj)
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
  * `i, j ::Int`
  * `Ri, Rj ::Vector{Int}`
"""
struct InteractionDiagonal{R<:Real}
    amplitude ::R
    i ::Int
    j ::Int
    Ri ::Vector{Int}
    Rj ::Vector{Int}
    function InteractionDiagonal{R}(v::R,
                                    i::Integer,
                                    j::Integer,
                                    Ri::AbstractVector{<:Integer},
                                    Rj::AbstractVector{<:Integer}) where {R<:Real}
        if length(Ri) != length(Rj)
            throw(ArgumentError("Ri and Rj must be of same length"))
        elseif (i, Ri) == (j, Rj)
            throw(ArgumentError("not allowed by fermion statistics"))
        end
        if (i, Ri) > (j, Rj)
            (i,j) = (j,i)
            (Ri, Rj) = (Rj, Ri)
            #v = conj(v)
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

Ordering of orbitals by `(i, Ri)`

  # Members
  * `amplitude ::C`
  * `i, j, k, l ::Int`
  * `Ri, Rj, Rk, Rl ::Vector{Int}`
"""
struct InteractionOffdiagonal{C<:Number}
    amplitude ::C
    i ::Int
    j ::Int
    k ::Int
    l ::Int
    Ri ::Vector{Int}
    Rj ::Vector{Int}
    Rk ::Vector{Int}
    Rl ::Vector{Int}
    function InteractionOffdiagonal{C}(v::C,
                                       i::Integer,
                                       j::Integer,
                                       k::Integer,
                                       l::Integer,
                                       Ri::AbstractVector{<:Integer},
                                       Rj::AbstractVector{<:Integer},
                                       Rk::AbstractVector{<:Integer},
                                       Rl::AbstractVector{<:Integer}) where {C<:Number}
        if !(length(Ri) == length(Rj) == length(Rk) == length(Rl))
            throw(ArgumentError("R's must be of same length"))
        end

        if (i, Ri) == (j, Rj)
            throw(ArgumentError("i, j not allowed by fermion statistics"))
        elseif (i, Ri) > (j, Rj)
            (i, j) = (j, i)
            (Ri, Rj) = (Rj, Ri)
            v = -v
        end

        if (k, Rk) == (l, Rl)
            throw(ArgumentError("k, l not allowed by fermion statistics"))
        elseif (k, Rk) > (l, Rl)
            (k, l) = (l, k)
            (Rk, Rl) = (Rl, Rk)
            v = -v
        end

        if (i, Ri, j, Rj) == (k, Rk, l, Rl)
            throw(ArgumentError("not an off-diagonal element"))
        elseif (i, Ri, j, Rj) > (k, Rk, l, Rl)
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
    hoppingbycarte{T}

Make a hopping element with cartesian coordinates.

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

Make a hopping element with cartesian coordinates.

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

Make an interaction element with cartesian coordinates.

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

Make an interaction element with cartesian coordinates.

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


"""
    islocal

Check if the hopping element is local (i.e. Ri is zero)

"""
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


"""
    localized

Return a hopping element that is local (i.e. Ri is zero)
"""
function localized(hopping::HoppingDiagonal{R}) where {R <: Real}
    return HoppingDiagonal{R}(hopping.amplitude, hopping.i, zero(hopping.Ri))
end

function localized(hopping::HoppingOffdiagonal{C}) where {C <: Number}
    return HoppingOffdiagonal{C}(hopping.amplitude,
                                 hopping.i,
                                 hopping.j,
                                 zero(hopping.Ri),
                                 hopping.Rj - hopping.Ri)
end

function localized(interaction::InteractionDiagonal{R}) where {R <: Real}
    return InteractionDiagonal{R}(interaction.amplitude,
                                  interaction.i,
                                  interaction.j,
                                  zero(interaction.Ri),
                                  interaction.Rj - interaction.Ri)
end

function localized(interaction::InteractionOffdiagonal{C}) where {C <: Number}
    # Since i,j,k,l already is ordered in the interaction term, not need to worry about which one to localize.
    return InteractionOffdiagonal{C}(interaction.amplitude,
                                     interaction.i,
                                     interaction.j,
                                     interaction.k,
                                     interaction.l,
                                     zero(interaction.Ri),
                                     interaction.Rj - interaction.Ri,
                                     interaction.Rk - interaction.Ri,
                                     interaction.Rl - interaction.Ri)
end

import Base.isapprox

function isapprox(lhs ::HoppingDiagonal{R},
                  rhs ::HoppingDiagonal{R};
                  atol ::Real=sqrt(eps(Float64)),
                  rtol ::Real=sqrt(eps(Float64))) where {R <: Real}
    return (lhs.i == rhs.i && lhs.Ri == rhs.Ri && isapprox(lhs.amplitude, rhs.amplitude; atol=atol, rtol=rtol))
end

function isapprox(lhs ::HoppingOffdiagonal{C},
                  rhs ::HoppingOffdiagonal{C};
                  atol ::Real=sqrt(eps(Float64)),
                  rtol ::Real=sqrt(eps(Float64))) where {C <:Number}
    return (lhs.i == rhs.i &&
            lhs.j == rhs.j &&
            lhs.Ri == rhs.Ri &&
            lhs.Rj == rhs.Rj &&
            isapprox(lhs.amplitude, rhs.amplitude; atol=atol, rtol=rtol))
end

function isapprox(lhs ::InteractionDiagonal{R},
                  rhs ::InteractionDiagonal{R};
                  atol ::Real=sqrt(eps(Float64)),
                  rtol ::Real=sqrt(eps(Float64))) where {R <: Real}
    return (lhs.i == rhs.i &&
            lhs.j == rhs.j &&
            lhs.Ri == rhs.Ri &&
            lhs.Rj == rhs.Rj &&
            isapprox(lhs.amplitude, rhs.amplitude; atol=atol, rtol=rtol))
end

function isapprox(lhs ::InteractionOffdiagonal{C},
                  rhs ::InteractionOffdiagonal{C};
                  atol ::Real=sqrt(eps(Float64)),
                  rtol ::Real=sqrt(eps(Float64))) where {C <: Number}
    return (lhs.i == rhs.i &&
            lhs.j == rhs.j &&
            lhs.k == rhs.k &&
            lhs.l == rhs.l &&
            lhs.Ri == rhs.Ri &&
            lhs.Rj == rhs.Rj &&
            lhs.Rk == rhs.Rk &&
            lhs.Rl == rhs.Rl &&
            isapprox(lhs.amplitude, rhs.amplitude; atol=atol, rtol=rtol))
end


import Base.isless

isless(hop1 ::HoppingDiagonal, hop2::HoppingOffdiagonal) = true
isless(hop1 ::HoppingOffdiagonal, hop2::HoppingDiagonal) = false
