__precompile__()
module HFB

using ..Lattice
using ..Spec
using ..Generator

using PyCall
#@pyimport numpy.linalg as npl

const npl = PyNULL()

function __init__()
    copy!(npl, pyimport_conda("numpy.linalg", "numpy"))
end

include("model.jl")
include("computer.jl")
include("solver.jl")
include("freeenergy.jl")

end
