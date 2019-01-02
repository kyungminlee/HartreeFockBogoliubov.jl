"""
    Submodule `Spec`

"""
module Spec

using ..Lattice

export HoppingDiagonal,
       HoppingOffdiagonal,
       InteractionDiagonal,
       InteractionOffdiagonal
export Hopping, Interaction
export hoppingbycarte, interactionbycarte
export addhopping!, addinteraction!
export islocal, localized

export FullHamiltonian

include("./terms.jl")
include("./hamiltonian.jl")

end # Spec
