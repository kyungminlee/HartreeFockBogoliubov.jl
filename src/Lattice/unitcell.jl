export UnitCell
export make_unitcell
export dimension,
       numorbital,
       hasorbital,
       addorbital!,
       getorbital,
       getorbitalindex,
       getorbitalcoord,
       getorbitalindexcoord,
       getorbitalname,
       carte2fract,
       fract2carte,
       whichunitcell,
       momentumgrid

using LinearAlgebra
import Base.zero


"""
    UnitCell{O}

# Parameters
* `O`: type of "orbital". Any type can be used, but we recommend using `String` or tuple of `String` and `Int`
       for compatibility with JSON.

# Members
* `latticevectors ::Array{Float64, 2}`: Lattice vectors
* `reducedreciprocallatticevectors ::Array{Float64, 2}`: Reduced reciprocal lattice vectors (transpose of inverse of `latticevectors`)
* `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors. `2π * reducedreciprocallatticevectors`
* `orbitals ::Vector{Tuple{T, FractCoord}}`: List of orbitals within unit cell
* `orbitalindices ::Dict{T, Int}`: Indices of orbitals
"""
mutable struct UnitCell{O}
    latticevectors ::Array{Float64, 2}
    orbitals ::Vector{Tuple{O, FractCoord}}

    reducedreciprocallatticevectors ::Array{Float64, 2}
    reciprocallatticevectors ::Array{Float64, 2}
    orbitalindices ::Dict{O, Int}
    function UnitCell{O}(latticevectors ::AbstractArray{<:Real, 2},
                         orbitals ::AbstractVector{Tuple{O, FractCoord}},
                         reducedreciprocallatticevectors ::AbstractArray{<:Real, 2},
                         reciprocallatticevectors ::AbstractArray{<:Real, 2},
                         orbitalindices ::Dict{O, Int}) where {O}
        if (O <: Integer)
            throw(ArgumentError("OrbitalType should not be integer to avoid confusion"))
        end
        new{O}(latticevectors,
               orbitals,
               reducedreciprocallatticevectors,
               reciprocallatticevectors,
               orbitalindices)
    end
end


"""
    UnitCell

Construct a one-dimensional lattice.

# Arguments
* `latticeconstant ::Float64`: Lattice constant
* `OrbitalType`: List of orbitals

# Optional Arguments
* `tol=sqrt(eps(Float64))`: Tolerance
"""
function make_unitcell(latticeconstant ::AbstractFloat;
                       OrbitalType::DataType=Any,
                       tol=sqrt(eps(Float64)))
    return make_unitcell(reshape([latticeconstant], (1,1));
                         OrbitalType=OrbitalType, tol=tol)
end


"""
    UnitCell

Construct an n-dimensional lattice.

# Arguments
* `latticevectors ::AbstractArray{<:AbstractFloat, 2}`: Lattice vectors
* `OrbitalType::DataType`

# Optional Arguments
* `tol=sqrt(eps(Float64))`: Epsilon
"""
function make_unitcell(latticevectors ::AbstractArray{<:AbstractFloat, 2};
                       OrbitalType::DataType=Any,
                       tol ::Real=sqrt(eps(Float64)))
    (ndim, ndim_) = size(latticevectors)
    if ndim != ndim_
        throw(ArgumentError("lattice vectors has dimension ($(ndim), $(ndim_))"))
    elseif abs(det(latticevectors)) <= tol
        throw(ArgumentError("lattice vectors define zero volume $(latticevectors)"))
    elseif tol < 0
        throw(ArgumentError("tol must be non-negative"))
    end

    reduced_rlv = transpose(inv(latticevectors))
    orbitals = Tuple{OrbitalType, FractCoord}[]
    orbitalindices = Dict{OrbitalType, Int}()
    return UnitCell{OrbitalType}(latticevectors, orbitals,
                                 reduced_rlv, 2*π*reduced_rlv, orbitalindices)
end


function make_unitcell(latticevectors::AbstractVector{<:AbstractVector};
                       OrbitalType::DataType=Any,
                       tol ::Real=sqrt(eps(Float64)))
    lv = hcat(latticevectors...)
    return make_unitcell(lv; OrbitalType=OrbitalType, tol=tol)
end


"""
    dimension

Spatial dimension of the unit cell.
"""
function dimension(uc ::UnitCell)
    return size(uc.latticevectors)[1]
end


"""
    numorbital

Number of orbitals of the unit cell.

# Arguments
* `uc ::UnitCell`
"""
function numorbital(uc ::UnitCell)
    return length(uc.orbitals)
end


"""
    addorbital!

Add an orbital to the unit cell.

# Arguments
* `uc ::UnitCell{T}`
* `orbitalname ::{T}`
* `orbitalcoord ::FractCoord`
"""
function addorbital!(uc ::UnitCell{O},
                     orbitalname ::O,
                     orbitalcoord ::FractCoord) where {O}
    (ndim, ndim_) = size(uc.latticevectors)
    if dimension(orbitalcoord) != ndim
        throw(ArgumentError("orbitalcoord has wrong dimension"))
    elseif haskey(uc.orbitalindices, orbitalname)
        throw(ArgumentError( "duplicate orbital name"))
    end
    push!(uc.orbitals, (orbitalname, orbitalcoord))
    index = length(uc.orbitalindices)+1
    uc.orbitalindices[orbitalname] = index
    return index
end


"""
    hasorbital{T}

Test whether the unit cell contains the orbital of given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`
"""
function hasorbital(uc ::UnitCell{O}, name ::O) where {O}
    return haskey(uc.orbitalindices, name)
end


"""
    getorbitalindex

Get index of the given orbital.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`
"""
function getorbitalindex(uc ::UnitCell{O}, name ::O) where {O}
    return uc.orbitalindices[name]
end


"""
    getorbital

    Get the orbital (its orbital name and its fractional coordinates) with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `(orbitalname, fractcoord)`
"""
function getorbital(uc ::UnitCell{O}, name ::O) where {O}
    return uc.orbitals[ uc.orbitalindices[name] ]
end


"""
    getorbitalcoord

    Get the fractional coordinates of the orbital with the given name.

# Arguments
* `uc ::UnitCell{O}`
* `name ::O`

# Return
* `fractcoord`
"""
function getorbitalcoord(uc ::UnitCell{O}, name ::O) where {O}
    return getorbital(uc, name)[2]
end


"""
    getorbitalindexcoord

# Arguments
* `uc ::UnitCell{T}`
* `name ::T`

# Return
* `(index, fractcoord)`
"""
function getorbitalindexcoord(uc ::UnitCell{O}, name::O) where {O}
    index = getorbitalindex(uc, name)
    coord = getorbitalcoord(uc, index)
    return (index, coord)
end


"""
    getorbital

# Arguments
* `uc ::UnitCell{T}`
* `index ::Integer`

# Return
* `(orbitalname, fractcoord)`
"""
function getorbital(uc ::UnitCell, index::Integer)
    return uc.orbitals[index]
end


"""
    getorbitalname

# Arguments
* `uc ::UnitCell`
* `index ::Integer`

# Return
* `orbitalname`
"""
function getorbitalname(uc ::UnitCell, index ::Integer)
    return uc.orbitals[index][1]
end


"""
    getorbitalcoord

# Arguments
* `uc ::UnitCell`
* `idx ::Integer`

# Return
* `fractcoord`
"""
function getorbitalcoord(uc ::UnitCell, index ::Integer)
    return uc.orbitals[index][2]
end



"""
    fract2carte

# Arguments
* `latticevectors ::Array{Float64, 2}`
* `fc ::FractCoord`
"""
function fract2carte(unitcell ::UnitCell, fc ::FractCoord)
    if dimension(unitcell) != dimension(fc)
        throw(ArgumentError("unitcell and fractcoord must have the same dimension"))
    end
    mc = fc.whole + fc.fraction
    cc = unitcell.latticevectors * mc
    return CarteCoord(cc)
end


"""
    carte2fract

# Arguments
* `latticevectors ::Array{Float64, 2}`
* `cc ::CarteCoord`
"""
function carte2fract(unitcell ::UnitCell,
                     cc ::CarteCoord;
                     tol ::Real=sqrt(eps(Float64)))
    if dimension(unitcell) != length(cc)
        throw(ArgumentError("unitcell and cartecoord must have the same dimension"))
    end
    fc = transpose(unitcell.reducedreciprocallatticevectors) * cc
    w = Int[fld(x, 1.0) for x in fc]
    f = Float64[mod(x, 1.0) for x in fc]
    for i in length(w)
        if f[i] + tol >= 1.0
            w[i] += 1
            f[i] = 0.0
        end
    end
    return FractCoord(w, f)
end


"""
    whichunitcell

# Return
* `R ::Vector{Int}`: which unit cell the specificied orbital/cartesian coordinates belongs to.
"""
function whichunitcell(uc ::UnitCell{O},
                       name ::O,
                       cc ::CarteCoord;
                       tol ::Real=sqrt(eps(Float64))) where {O}
    fc1 = getorbitalcoord(uc, name)
    fc2 = carte2fract(uc, cc)
    if !isapprox(fc1.fraction, fc2.fraction; rtol=tol)
        throw(ArgumentError("$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]"))
    end
    R = fc2.whole - fc1.whole
    return R
end


"""
    whichunitcell

# Return
* `R ::Vector{Int}`: which unit cell the specificied orbital/cartesian coordinates belongs to.
"""
function whichunitcell(uc ::UnitCell{O},
                       name ::O,
                       fc ::FractCoord;
                       tol ::Real=sqrt(eps(Float64))) where {O}
    fc1 = getorbitalcoord(uc, name)
    fc2 = fc
    if !isapprox(fc1.fraction, fc2.fraction; rtol=tol)
        throw(ArgumentError("$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]"))
    end
    R = fc2.whole - fc1.whole
    return R
end


function zero(uc::UnitCell; dtype::DataType=ComplexF64)
    norb = numorbital(uc)
    return Base.zeros(dtype, (norb, norb))
end


"""
    momentumgrid

Generate an n-dimensional grid of momenta of given shape
"""
function momentumgrid(uc::UnitCell, shape::AbstractVector{<:Integer})
    if length(shape) != dimension(uc)
        throw(ArgumentError("dimension mismatch"))
    elseif !all((x) -> x>0, shape)
        throw(ArgumentError("shape should be positive"))
    end
    ranges = [range(0,stop=1,length=n+1)[1:end-1] for n in shape]
    cubicgrid = map((x) -> [x...], Base.product(ranges...))
    momentumgrid = map((x) -> uc.reciprocallatticevectors * x, cubicgrid)
    return momentumgrid
end
