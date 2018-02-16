__precompile__(false)  ## because HFB uses pyimport
module HFB

using ..Lattice
using ..Spec
using ..Generator

include("model.jl")
include("computer.jl")
include("solver.jl")
include("freeenergy.jl")

end
