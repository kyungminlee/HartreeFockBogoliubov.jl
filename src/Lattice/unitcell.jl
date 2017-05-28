
export UnitCell

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

  """
      UnitCell

    # Arguments
    * `latticevectors ::Array{Float64, 2}`: Lattice vectors
    * `orbitals::Vector{Tuple{AbstractString, FractCoord}}=[]`: List of orbitals

    # Optional Arguments
    * `epsilon=sqrt(eps(Float64))`: Epsilon
  """
  function UnitCell(latticevectors ::Array{Float64, 2}, 
                    orbitals::Vector{Tuple{AbstractString, FractCoord}}=[];
                    epsilon=sqrt(eps(Float64)))
    (ndim, ndim_) = size(latticevectors)
    @assert(ndim -= ndim_, "lattice vectors has dimension ($(ndim), $(ndim_))")
    @assert(det(latticevectors) > epsilon, "lattice vectors define zero volume")

    reduced_rlv = inv(transpose(latticevectors))

    if isempty(orbitals)
      orbitals = [("A", FractCoord(ndim))]
      orbitalindices = Dict{AbstractString, Int64}("A" => 1)
      return new(latticevectors, reduced_rlv, 2*pi*reduced_rlv, orbitals, orbitalindices)
    else
      orbitalindices = Dict{AbstractString, Int64}()
      for (idx, orbital) in enumerate(orbitals)
        @assert(!haskey(orbitalindices, orbital[1]), "UnitCell: duplicate orbital name $(orbital[1])")
        @assert(dimension(orbital[2]) == ndim, "orbital coordinate has wrong dimension")
        orbitalindices[orbital[1]] = idx
      end
      return new(latticevectors, reduced_rlv, 2*pi*reduced_rlv, orbitals, orbitalindices)
    end
  end
end


"""
    addorbital!

  # Arguments
  * `uc ::UnitCell`
  * `orbitalname ::AbstractString`
  * `orbitalcoord ::FractCoord`
"""
function addorbital!(uc ::UnitCell, orbitalname ::AbstractString, orbitalcoord ::FractCoord)
  (ndim, ndim_) = size(uc.latticevectors)
  @assert(uc.latticevectors[1] == dimension(orbitalcoord), "orbitalcoord has wrong dimension")
  @assert(!haskey(uc.orbital_to_coord, orbitalname), "duplicate orbital name")
  push!(uc.orbitals, (orbitalname, orbitalcoord))
  uc.orbitalindices[orbitalname] = length(uc.orbitalindices)
end

"""
    hasorbital

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function hasorbital(uc ::UnitCell, name ::AbstractString)
  return haskey(uc.orbitals, name)
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
    getorbitalindex

  # Arguments
  * `uc ::UnitCell`
  * `name ::AbstractString`
"""
function getorbitalindex(uc ::UnitCell, name ::AbstractString)
  return uc.orbitalindices[name]
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
    numorbital

  # Arguments
  * `uc ::UnitCell`
"""
function numorbital(uc ::UnitCell)
  return length(uc.orbitals)
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


