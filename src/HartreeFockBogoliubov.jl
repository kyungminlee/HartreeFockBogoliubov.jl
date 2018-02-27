__precompile__()
module HartreeFockBogoliubov
include("basic.jl")

include("Lattice/Lattice.jl")
using .Lattice

for symb in names(Lattice)
    eval(Expr(:export, symb))
end

include("Hamiltonian/Spec.jl")
include("Hamiltonian/Generator.jl")

include("HFB/HFB.jl")

include("Topology/Topology.jl")

include("LinearizedGap/LinearizedGap.jl")

include("IO/dumper.jl")
include("IO/dictify.jl")
include("IO/loader.jl")

end
