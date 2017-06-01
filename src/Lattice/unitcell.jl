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
       whichunitcell


"""
    UnitCell{T}

  # Members
  * `latticevectors ::Array{Float64, 2}`: Lattice vectors
  * `reducedreciprocallatticevectors ::Array{Float64, 2}`: Reduced reciprocal lattice vectors (transpose of inverse of `latticevectors`)
  * `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors
  * `orbitals ::Vector{Tuple{T, FractCoord}}`: List of orbitals within unit cell
  * `orbitalindices ::Dict{T, Int64}`: Indices of orbitals
"""
type UnitCell{T}
  latticevectors ::Array{Float64, 2}
  reducedreciprocallatticevectors ::Array{Float64, 2}
  reciprocallatticevectors ::Array{Float64, 2}

  orbitals ::Vector{Tuple{T, FractCoord}}
  orbitalindices ::Dict{T, Int64}
end


"""
    UnitCell

  Construct a one-dimensional lattice.

  # Arguments
  * `latticeconstant ::Float64`: Lattice constant
  * `OrbitalType`: List of orbitals

  # Optional Arguments
  * `tol=sqrt(eps(Float64))`: Epsilon
"""
function newunitcell(latticeconstant ::Float64;
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
function newunitcell(latticevectors ::Array{Float64, 2};
                     OrbitalType::DataType=Any,
                     tol=sqrt(eps(Float64)))
  (ndim, ndim_) = size(latticevectors)
  @assert(ndim == ndim_, "lattice vectors has dimension ($(ndim), $(ndim_))")
  @assert(abs(det(latticevectors)) > tol,
          "lattice vectors define zero volume $(latticevectors)")

  reduced_rlv = transpose(inv(latticevectors))
  orbitals = Tuple{OrbitalType, FractCoord}[]
  orbitalindices = Dict{OrbitalType, Int64}()
  return UnitCell{OrbitalType}(latticevectors, reduced_rlv, 2*pi*reduced_rlv, orbitals, orbitalindices)
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
function addorbital!{T}(uc ::UnitCell{T}, orbitalname ::T, orbitalcoord ::FractCoord)
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
  * `uc ::UnitCell{T}`
  * `name ::T`
"""
function hasorbital{T}(uc ::UnitCell{T}, name ::T)
  return haskey(uc.orbitalindices, name)
end

"""
    getorbitalindex

  Get index of the given orbital.

  # Arguments
  * `uc ::UnitCell{T}`
  * `name ::T`
"""
function getorbitalindex{T}(uc ::UnitCell{T}, name ::T)
  return uc.orbitalindices[name]
end

"""
    getorbital

  Get the orbital with the given name.

  # Arguments
  * `uc ::UnitCell{T}`
  * `name ::T`
"""
function getorbital{T}(uc ::UnitCell{T}, name ::T)
  return uc.orbitals[ uc.orbitalindices[name] ]
end

"""
    getorbitalcoord

  # Arguments
  * `uc ::UnitCell{T}`
  * `name ::T`
"""
function getorbitalcoord{T}(uc ::UnitCell{T}, name ::T)
  return getorbital(uc, name)[2]
end

"""
    getorbitalindexcoord

  # Arguments
  * `uc ::UnitCell{T}`
  * `name ::T`
"""
function getorbitalindexcoord{T}(uc ::UnitCell{T}, name::T)
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
function carte2fract(unitcell ::UnitCell, cc ::CarteCoord)
  fc = inv(unitcell.latticevectors) * cc
  w = Int64[fld(x, 1) for x in fc]
  f = Float64[mod(x, 1) for x in fc]
  return FractCoord(w, f)
end



"""
    whichunitcell

  Return which unit cell the specificied orbital/cartesian coordinates belongs to.
"""
function whichunitcell{T}(uc ::UnitCell{T}, name ::T, cc ::CarteCoord; tol=sqrt(eps(Float64)))
  fc1 = getorbitalcoord(uc, name)
  fc2 = carte2fract(uc, cc)
  @assert isapprox(fc1.fraction, fc2.fraction; rtol=tol)
  R = fc2.whole - fc1.whole
  return R
end
