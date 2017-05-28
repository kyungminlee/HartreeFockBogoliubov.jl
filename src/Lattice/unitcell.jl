
export UnitCell

type UnitCell
  latticevectors ::Array{Float64, 2}
  reducedreciprocallatticevectors ::Array{Float64, 2}
  reciprocallatticevectors ::Array{Float64, 2}

  orbitals ::Vector{Tuple{AbstractString, FractCoord}}
  orbitalindices ::Dict{AbstractString, Int64}

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


function addorbital!(uc ::UnitCell, orbitalname ::AbstractString, orbitalcoord ::FractCoord)
  (ndim, ndim_) = size(uc.latticevectors)
  @assert(uc.latticevectors[1] == dimension(orbitalcoord), "orbitalcoord has wrong dimension")
  @assert(!haskey(uc.orbital_to_coord, orbitalname), "duplicate orbital name")
  push!(uc.orbitals, (orbitalname, orbitalcoord))
  uc.orbitalindices[orbitalname] = length(uc.orbitalindices)
end


function hasorbital(uc ::UnitCell, name ::AbstractString)
  return haskey(uc.orbitals, name)
end


function getorbital(uc ::UnitCell, name ::AbstractString)
  return uc.orbitals[ uc.orbitalindices[name] ]
end


function getorbitalindex(uc ::UnitCell, name ::AbstractString)
  return uc.orbitalindices[name]
end


function getorbitalcoord(uc ::UnitCell, name ::AbstractString)
  return getorbital(uc, name)[2]
end


function getorbitalindexcoord(uc ::UnitCell, name::AbstractString)
  idx = getorbitalindex(uc, name)
  coord = getorbitalcoord(uc, idx)
  return (idx, coord)
end


function numorbital(uc ::UnitCell)
  return length(uc.orbitals)
end


function getorbital(uc ::UnitCell, idx::Integer)
  return uc.orbitals[idx]
end


function getorbitalname(uc ::UnitCell, idx ::Integer)
  return uc.orbitals[idx][1]
end


function getorbitalcoord(uc ::UnitCell, idx ::Integer)
  return uc.orbitals[idx][2]
end


function fract_to_carte(uc ::UnitCell, fc ::FractCoord)
  mc = fc.whole + fc.fraction
  cc = uc.latticevectors * mc
  return CarteCoord(cc)
end


function carte_to_fract(uc ::UnitCell, cc ::CarteCoord)
  fc = transpose(uc.reducedreciprocallatticevectors) * cc
  w = Int64[fld(x, 1) for x in fc]
  f = Float64[mod(x, 1) for x in fc]
  return FractCoord(w, f)
end

