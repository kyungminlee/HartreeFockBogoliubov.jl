export CarteCoord
export FractCoord
export dimension

"""
    CarteCoord

  Cartesian coordinates
"""
typealias CarteCoord Vector{Float64}


"""
    FractCoord

  Fractional coordinates
"""
type FractCoord
  whole ::Vector{Int64}
  fraction ::Vector{Float64}
  
  function FractCoord(coord ::Vector{Float64})
    w = Int64[fld(x,1) for x in coord]
    f = Float64[mod(x,1) for x in coord]
    return new(w, f)
  end

  function FractCoord(ndim ::Integer)
    @assert(ndim > 0, "ndim should be positive")
    w = zeros(Int64, ndim)
    f = zeros(Float64, ndim)
    return new(w, f)
  end
end

"""
    dimension

  Dimension of the fractional coordinates

  # Arguments
  * `fc ::FractCoord`: Fractional coordinates.
"""
function dimension(fc ::FractCoord)
  d1 = length(fc.whole)
  d2 = length(fc.fraction)
  @assert(d1 == d2, "FractCoord: wrong dimension")
  return d1
end


function +(fractcoord ::FractCoord, R ::Vector{Int64})
  return FractCoord(fractcoord.whole + R, fractcoord[2])
end


function -(fractcoord ::FractCoord, R ::Vector{Int64})
  return FractCoord(fractcoord.whole - R, fractcoord[2])
end


function +(fc1 ::FractCoord, fc2 ::FractCoord)
  R = [fld(x+y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  r = [mod(x+y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  return FractCoord(fc1.whole + fc2.whole + R, r)
end


function -(fc1 ::FractCoord, fc2 ::FractCoord)
  R = [fld(x-y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  r = [mod(x-y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  return FractCoord(fc1.whole - fc2.whole + R, r)
end



"""
    fract_to_carte
"""
function fract_to_carte(latticevectors ::Array{Float64, 2}, fc ::FractCoord)
  mc = fc.whole + fc.fraction
  cc = latticevectors * mc
  return CarteCoord(cc)
end


"""
    carte_to_fract
"""
function carte_to_fract(latticevectors ::Array{Float64, 2}, cc ::CarteCoord)
  fc = inv(latticevectors) * cc
  w = Int64[fld(x, 1) for x in fc]
  f = Float64[mod(x, 1) for x in fc]
  return FractCoord(w, f)
end

