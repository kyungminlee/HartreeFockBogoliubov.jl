## Basic functions

export squarify

"""
    squarify

In order to make the Hamiltoinian a periodic function of momentum, bring all the sites to the origin. In addition, make the unitcell into a square.

# Arguments
* `uc::Lattice.UnitCell{O}`
"""
function squarify(uc::Lattice.UnitCell{O}) where {O}
    newdimension = dimension(uc)
    newlatticevectors = eye(newdimension)

    origin = FractCoord(zeros(Int, newdimension), zeros(Float64, newdimension))

    newuc = make_unitcell(newlatticevectors; OrbitalType=O)
    for (orbname, fc) in uc.orbitals
        addorbital!(newuc, orbname, origin)
    end
    return newuc
end


"""
    squarify

# Arguments
* `uc::Spec.FullHamiltonian{O}`
"""
function squarify(ham::Spec.FullHamiltonian{O}) where {O}
    newuc = squarify(ham.unitcell)
    return FullHamiltonian{O}(newuc, ham.hoppings, ham.interactions)
end


"""
    squarify

# Arguments
* `uc::HFB.HFBHamiltonian{O}`
"""
function squarify(ham::HFB.HFBHamiltonian{O}) where {O}
    newuc = squarify(ham.unitcell)
    return HFBHamiltonian{O}(newuc,
                             ham.hoppings,
                             ham.particle_hole_interactions,
                             ham.particle_particle_interactions)
end
