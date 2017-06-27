export UnitCell
export newunitcell
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

import Base.zeros


"""
    UnitCell{T}

  # Members
  * `latticevectors ::Array{Float64, 2}`: Lattice vectors
  * `reducedreciprocallatticevectors ::Array{Float64, 2}`: Reduced reciprocal lattice vectors (transpose of inverse of `latticevectors`)
  * `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors
  * `orbitals ::Vector{Tuple{T, FractCoord}}`: List of orbitals within unit cell
  * `orbitalindices ::Dict{T, Int64}`: Indices of orbitals
"""
mutable struct UnitCell{O}
  latticevectors ::Array{Float64, 2}
  orbitals ::Vector{Tuple{O, FractCoord}}

  reducedreciprocallatticevectors ::Array{Float64, 2}
  reciprocallatticevectors ::Array{Float64, 2}
  orbitalindices ::Dict{O, Int64}
  function UnitCell{O}(latticevectors ::AbstractArray{<:Real, 2},
                       orbitals ::AbstractVector{Tuple{O, FractCoord}},
                       reducedreciprocallatticevectors ::AbstractArray{<:Real, 2},
                       reciprocallatticevectors ::AbstractArray{<:Real, 2},
                       orbitalindices ::Dict{O, Int}) where {O}
    @assert(! (O <: Integer), "OrbitalType should not be integer to avoid confusion")
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
function newunitcell(latticeconstant ::AbstractFloat;
                     OrbitalType::DataType=Any,
                     tol=sqrt(eps(Float64)))
  return newunitcell(reshape([latticeconstant], (1,1));
                     OrbitalType=OrbitalType, tol=tol)
end


"""
    UnitCell

  # Arguments
  * `latticevectors ::Array{Float64, 2}`: Lattice vectors
  * `OrbitalType::DataType`

  # Optional Arguments
  * `tol=sqrt(eps(Float64))`: Epsilon
"""
function newunitcell(latticevectors ::AbstractArray{<:AbstractFloat, 2};
                     OrbitalType::DataType=Any,
                     tol=sqrt(eps(Float64)))
  (ndim, ndim_) = size(latticevectors)
  @assert(ndim == ndim_, "lattice vectors has dimension ($(ndim), $(ndim_))")
  @assert(abs(det(latticevectors)) > tol,
          "lattice vectors define zero volume $(latticevectors)")

  reduced_rlv = transpose(inv(latticevectors))
  orbitals = Tuple{OrbitalType, FractCoord}[]
  orbitalindices = Dict{OrbitalType, Int64}()
  return UnitCell{OrbitalType}(latticevectors, orbitals,
                               reduced_rlv, 2*pi*reduced_rlv, orbitalindices)
end


function newunitcell(latticevectors::AbstractVector{<:AbstractVector};
                     OrbitalType::DataType=Any,
                     tol=sqrt(eps(Float64)))
  lv = hcat(latticevectors...)
  return newunitcell(lv; OrbitalType=OrbitalType, tol=tol)
end

"""
    dimension

  Spatial dimension of the unit cell.
"""
function dimension(uc ::UnitCell{O}) where {O}
  return size(uc.latticevectors)[1]
end


"""
    numorbital

  Number of orbitals of the unit cell.

  # Arguments
  * `uc ::UnitCell`
"""
function numorbital(uc ::UnitCell{O}) where {O}
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
                     orbitalname ::O, orbitalcoord ::FractCoord) where {O}
  (ndim, ndim_) = size(uc.latticevectors)
  @assert(dimension(orbitalcoord) == ndim, "orbitalcoord has wrong dimension")
  @assert(!haskey(uc.orbitalindices, orbitalname), "duplicate orbital name")
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

  Get the orbital with the given name.

  # Arguments
  * `uc ::UnitCell{O}`
  * `name ::O`
"""
function getorbital(uc ::UnitCell{O}, name ::O) where {O}
  return uc.orbitals[ uc.orbitalindices[name] ]
end


"""
    getorbitalcoord

  # Arguments
  * `uc ::UnitCell{O}`
  * `name ::O`
"""
function getorbitalcoord(uc ::UnitCell{O}, name ::O) where {O}
  return getorbital(uc, name)[2]
end


"""
    getorbitalindexcoord

  # Arguments
  * `uc ::UnitCell{T}`
  * `name ::T`
"""
function getorbitalindexcoord(uc ::UnitCell{O}, name::O) where {O}
  idx = getorbitalindex(uc, name)
  coord = getorbitalcoord(uc, idx)
  return (idx, coord)
end


"""
    getorbital

  # Arguments
  * `uc ::UnitCell{T}`
  * `idx ::Integer`
"""
function getorbital(uc ::UnitCell, idx::Integer)
  return uc.orbitals[idx]
end


"""
    getorbitalname

  # Arguments
  * `uc ::UnitCell`
  * `idx ::Integer`
"""
function getorbitalname(uc ::UnitCell, idx ::Integer)
  return uc.orbitals[idx][1]
end


"""
    getorbitalcoord

  # Arguments
  * `uc ::UnitCell`
  * `idx ::Integer`
"""
function getorbitalcoord(uc ::UnitCell, idx ::Integer)
  return uc.orbitals[idx][2]
end



"""
    fract2carte

  # Arguments
  * `latticevectors ::Array{Float64, 2}`
  * `fc ::FractCoord`
"""
function fract2carte(unitcell ::UnitCell, fc ::FractCoord)
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
                     tol::Real=sqrt(eps(Float64)))
  #fc = inv(unitcell.latticevectors) * cc
  fc = transpose(unitcell.reducedreciprocallatticevectors) * cc
  w = Int64[fld(x, 1.0) for x in fc]
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

  Return which unit cell the specificied orbital/cartesian coordinates belongs to.
"""
function whichunitcell(uc ::UnitCell{O}, name ::O, cc ::CarteCoord;
                       tol::Real=sqrt(eps(Float64))) where {O}
  fc1 = getorbitalcoord(uc, name)
  fc2 = carte2fract(uc, cc)
  @assert(isapprox(fc1.fraction, fc2.fraction; rtol=tol),
          "$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]")
  R = fc2.whole - fc1.whole
  return R
end


function whichunitcell(uc ::UnitCell{O}, name ::O, fc ::FractCoord;
                       tol::Real=sqrt(eps(Float64))) where {O}
  fc1 = getorbitalcoord(uc, name)
  fc2 = fc
  @assert(isapprox(fc1.fraction, fc2.fraction; rtol=tol),
          "$(fc1.fraction) != $(fc2.fraction) [tol=$(tol)]")
  R = fc2.whole - fc1.whole
  return R
end


function zeros(uc::UnitCell; dtype::DataType=Complex128)
  norb = numorbital(uc)
  Base.zeros(dtype, (norb, norb))
end


function momentumgrid(uc::UnitCell, shape::AbstractVector{<:Integer})
  @assert(length(shape) == dimension(uc), "dimension mismatch")
  @assert(all((x) -> x>0, shape), "shape should be positive")

  ranges = [linspace(0,1,n+1)[1:end-1] for n in shape]

  cubicgrid = map((x) -> [x...], Base.product(ranges...))
  momentumgrid = map((x) -> transpose(uc.reciprocallatticevectors) * x, cubicgrid)
  return momentumgrid
end
