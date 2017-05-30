export UnitCell
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
    UnitCell

  # Members
  * `latticevectors ::Array{Float64, 2}`: Lattice vectors
  * `reducedreciprocallatticevectors ::Array{Float64, 2}`: Reduced reciprocal lattice vectors (transpose of inverse of `latticevectors`)
  * `reciprocallatticevectors ::Array{Float64, 2}`: Reciprocal lattice vectors
  * `orbitals ::Vector{Tuple{AbstractString, FractCoord}}`: List of orbitals within unit cell
  * `orbitalindices ::Dict{AbstractString, Int64}`: Indices of orbitals
"""
type UnitCell
  latticevectors ::Array{Float64, 2}
  reducedreciprocallatticevectors ::Array{Float64, 2}
  reciprocallatticevectors ::Array{Float64, 2}

  orbitals ::Vector{Tuple{AbstractString, FractCoord}}
  orbitalindices ::Dict{AbstractString, Int64}
end


"""
    UnitCell

  # Arguments
  * `latticevectors ::Array{Float64, 1}`: Lattice vectors
  * `orbitals::Vector{Tuple{AbstractString, FractCoord}}=[]`: List of orbitals

  # Optional Arguments
  * `epsilon=sqrt(eps(Float64))`: Epsilon
"""
function UnitCell(latticevectors ::Array{Float64, 1}, 
                  orbitals::Vector{Tuple{AbstractString, FractCoord}}
                    =Vector{Tuple{AbstractString, FractCoord}}();
                  tol=sqrt(eps(Float64)))
  @assert(size(latticevectors) == (1,), "lattice vectors should be 1D")
  return UnitCell(reshape(latticevectors, (1,1)), orbitals; tol=tol)
end


"""
    UnitCell

  # Arguments
  * `latticevectors ::Array{Float64, 2}`: Lattice vectors
  * `orbitals::Vector{Tuple{AbstractString, FractCoord}}=[]`: List of orbitals

  # Optional Arguments
  * `epsilon=sqrt(eps(Float64))`: Epsilon
"""
function UnitCell(latticevectors ::Array{Float64, 2}, 
                  orbitals::Vector{Tuple{AbstractString, FractCoord}}
                    =Vector{Tuple{AbstractString, FractCoord}}();
                  tol=sqrt(eps(Float64)))
  (ndim, ndim_) = size(latticevectors)
  @assert(ndim == ndim_, "lattice vectors has dimension ($(ndim), $(ndim_))")
  @assert(abs(det(latticevectors)) > tol, "lattice vectors define zero volume $(latticevectors)")

  reduced_rlv = transpose(inv(latticevectors))

  orbitalindices = Dict{AbstractString, Int64}()
  for (idx, orbital) in enumerate(orbitals)
    @assert(!haskey(orbitalindices, orbital[1]), "UnitCell: duplicate orbital name $(orbital[1])")
    @assert(dimension(orbital[2]) == ndim, "orbital coordinate has wrong dimension")
    orbitalindices[orbital[1]] = idx
  end

  return UnitCell(latticevectors, reduced_rlv, 2*pi*reduced_rlv, orbitals, orbitalindices)
end

"""
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

  Add orbital to the unit cell.

  # Arguments
  * `uc ::UnitCell`
  * `orbitalname ::AbstractString`
  * `orbitalcoord ::FractCoord`
"""
function addorbital!(uc ::UnitCell, orbitalname ::AbstractString, orbitalcoord ::FractCoord)
  (ndim, ndim_) = size(uc.latticevectors)
  @assert(dimension(orbitalcoord) == ndim, "orbitalcoord has wrong dimension")
  @assert(!haskey(uc.orbitalindices, orbitalname), "duplicate orbital name")
  push!(uc.orbitals, (orbitalname, orbitalcoord))
  index = length(uc.orbitalindices)+1
  uc.orbitalindices[orbitalname] = index
  return index
end

"""
    hasorbital

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function hasorbital(uc ::UnitCell, name ::AbstractString)
  return haskey(uc.orbitalindices, name)
end

"""
    getorbitalindex

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function getorbitalindex(uc ::UnitCell, name ::AbstractString)
  return uc.orbitalindices[name]
end

"""
    getorbital

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function getorbital(uc ::UnitCell, name ::AbstractString)
  return uc.orbitals[ uc.orbitalindices[name] ]
end

"""
    getorbitalcoord

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function getorbitalcoord(uc ::UnitCell, name ::AbstractString)
  return getorbital(uc, name)[2]
end

"""
    getorbitalindexcoord

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function getorbitalindexcoord(uc ::UnitCell, name::AbstractString)
  idx = getorbitalindex(uc, name)
  coord = getorbitalcoord(uc, idx)
  return (idx, coord)
end


"""
    getorbital

  # Arguments
  * `uc ::UnitCell`
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
function whichunitcell(uc ::UnitCell, name ::AbstractString, cc ::CarteCoord)
  fc1 = getorbitalcoord(uc, name)
  fc2 = carte2fract(uc, cc)
  @assert isapprox(fc1.fraction, fc2.fraction)
  R = fc2.whole - fc1.whole
  return R
end