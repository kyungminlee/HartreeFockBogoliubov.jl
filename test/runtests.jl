using Test

using HartreeFockBogoliubov
using HartreeFockBogoliubov.Spec
using HartreeFockBogoliubov.Generator
using HartreeFockBogoliubov.HFB


@generated function ≂(x, y)
    if !isempty(fieldnames(x)) && x == y
        mapreduce(n -> :(x.$n == y.$n), (a,b)->:($a && $b), fieldnames(x))
    else
        :(x == y)
    end
end

include("coord.jl")
include("unitcell.jl")

include("spec.jl")
include("generator.jl")

#include("hfb.jl")
#include("dictify.jl")
