export CarteCoord
export FractCoord
export dimension
export fract2carte, carte2fract


"""
    CarteCoord

  Cartesian coordinates. `Vector{Float64}`.
"""
const CarteCoord = Vector{Float64}


"""
    FractCoord

  Fractional coordinates.

  # Members
  * `whole ::Vector{Int64}`: Integer part of fractional coordinates
  * `fraction ::Vector{Float64}`: [0,1) part of fractional coordinates
"""
immutable FractCoord
  whole ::Vector{Int64}
  fraction ::Vector{Float64}

  function FractCoord(coord ::Vector{Float64})
    w = Int64[fld(x,1) for x in coord]
    f = Float64[mod(x,1) for x in coord]
    return new(w, f)
  end

  function FractCoord(w ::Vector{Int64}, f ::Vector{Float64})
    # TODO(kmlee) check
    @assert length(w) == length(f)
    @assert all(x -> 0 <= x < 1, f)
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

import Base: +, -


function -(fc ::FractCoord)
  return FractCoord(-(fc.whole + fc.fraction))
end


function +(fractcoord ::FractCoord, R ::Vector{Int64})
  return FractCoord(fractcoord.whole + R, fractcoord.fraction)
end


function -(fractcoord ::FractCoord, R ::Vector{Int64})
  return FractCoord(fractcoord.whole - R, fractcoord.fraction)
end


function +(fc1 ::FractCoord, fc2 ::FractCoord)
  R = Int64[fld(x+y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  r = Float64[mod(x+y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  return FractCoord(fc1.whole + fc2.whole + R, r)
end


function -(fc1 ::FractCoord, fc2 ::FractCoord)
  R = Int64[fld(x-y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  r = Float64[mod(x-y,1) for (x,y) in zip(fc1.fraction, fc2.fraction)]
  return FractCoord(fc1.whole - fc2.whole + R, r)
end

import Base.isapprox

function isapprox(fc1 ::FractCoord, fc2 ::FractCoord;
  rtol::Float64=sqrt(eps(Float64)) )
  return (fc1.whole == fc2.whole) &&
         isapprox(fc1.fraction, fc2.fraction; rtol=rtol)
end

import Base.show
function show(io::IO, fc::FractCoord)
  print(io, "FractCoord(", fc.whole, " + ", fc.fraction,")")
end



"""
    fract2carte

  # Arguments
  * `latticevectors ::Array{Float64, 2}`
  * `fc ::FractCoord`
"""
function fract2carte(latticevectors ::Array{Float64, 2}, fc ::FractCoord)
  mc = fc.whole + fc.fraction
  cc = latticevectors * mc
  return CarteCoord(cc)
end


"""
    carte2fract

  # Arguments
  * `latticevectors ::Array{Float64, 2}`
  * `cc ::CarteCoord`
"""
function carte2fract(latticevectors ::Array{Float64, 2}, cc ::CarteCoord)
  fc = inv(latticevectors) * cc
  w = Int64[fld(x, 1) for x in fc]
  f = Float64[mod(x, 1) for x in fc]
  return FractCoord(w, f)
end
