__precompile__()
module Dictify

using DataStructures

using ..Lattice
import ..Spec
import ..HFB

export dictify, objectify

dictify(obj::Number) = obj
dictify(obj::AbstractString) = obj
dictify(obj::AbstractArray) = [dictify(x) for x in obj]
dictify(obj::Dict) = Dict(dictify(k) => dictify(v) for (k, v) in obj)
dictify(obj::Symbol) = OrderedDict("type" => "Symbol", "value" => string(obj))

function dictify(obj::Tuple)
    OrderedDict(
        "type" => "Tuple",
        "value" => [dictify(x) for x in obj])
end

function dictify(fc::FractCoord)
    OrderedDict(
        "type" => "FractCoord",
        "whole" => fc.whole,
        "fraction" => fc.fraction)
end

function dictify(uc::UnitCell{O}) where {O}
    OrderedDict(
        "type" => "UnitCell",
        "latticevectors" => uc.latticevectors,
        "orbitals" => dictify(uc.orbitals))
end

function dictify(hopdiag::Spec.HoppingDiagonal{R}) where {R}
    OrderedDict(
        "type" => "HoppingDiagonal",
        "amplitude" => hopdiag.amplitude,
        "i" => hopdiag.i,
        "Ri" => hopdiag.Ri)
end

function dictify(hopoff::Spec.HoppingOffdiagonal{C}) where {C}
    OrderedDict(
        "type" => "HoppingOffdiagonal",
        "amplitude" => hopoff.amplitude,
        "i" => hopoff.i,
        "j" => hopoff.j,
        "Ri" => hopoff.Ri,
        "Rj" => hopoff.Rj,
        )
end

function dictify(intdiag ::Spec.InteractionDiagonal{R}) where {R}
    OrderedDict(
        "type" => "InteractionDiagonal",
        "amplitude" => intdiag.amplitude,
        "i" => intdiag.i,
        "j" => intdiag.j,
        "Ri" => intdiag.Ri,
        "Rj" => intdiag.Rj,
        )
end

function dictify(intoff ::Spec.InteractionOffdiagonal{C}) where {C}
    OrderedDict(
        "type" => "InteractionOffdiagonal",
        "amplitude" => intoff.amplitude,
        "i" => intoff.i,
        "j" => intoff.j,
        "k" => intoff.k,
        "l" => intoff.l,
        "Ri" => intoff.Ri,
        "Rj" => intoff.Rj,
        "Rk" => intoff.Rk,
        "Rl" => intoff.Rl,
        )
end

function dictify(hamspec::Spec.FullHamiltonian{O}) where {O}
    OrderedDict(
        "type" => "FullHamiltonian",
        "unitcell" => dictify(hamspec.unitcell),
        "hoppings" => dictify(hamspec.hoppings),
        "interactions" => dictify(hamspec.interactions)
        )
end

function dictify(hop::HFB.HoppingMeanField{R}) where {R}
    OrderedDict(
        "type" => "HoppingMeanField",
        "amplitude" => hop.amplitude,
        "target" => dictify(hop.target),
        "source" => dictify(hop.source),
        "star" => hop.star,
        )
end

function dictify(pairing::HFB.PairingMeanField{R}) where {R}
    OrderedDict(
        "type" => "PairingMeanField",
        "amplitude" => pairing.amplitude,
        "target" => dictify(pairing.target),
        "source" => dictify(pairing.source),
        "negate" => pairing.negate,
        )
end

function dictify(ham::HFB.HFBHamiltonian{O}) where {O}
    OrderedDict(
        "type" => "HFBHamiltonian",
        "unitcell" => dictify(ham.unitcell),
        "hoppings" => dictify(ham.hoppings),
        "particle_hole_interactions" => dictify(ham.particle_hole_interactions),
        "particle_particle_interactions" => dictify(ham.particle_particle_interactions),
        )
end

function dictify(computer::HFB.HFBComputer{O}) where {O}
    OrderedDict(
        "type" => "HFBComputer",
        "unitcell" => dictify(computer.unitcell),
        "hoppings" => dictify(computer.hoppings),
        "temperature" => computer.temperature,
        "rho_registry" => dictify(computer.ρ_registry),
        "t_registry"   => dictify(computer.t_registry),
        "Gamma_registry" => dictify(computer.Γ_registry),
        "Delta_registry" => dictify(computer.Δ_registry),
    )
end


function dictify(solver::HFB.HFBSolver{O}) where {O}
    OrderedDict(
        "type" => "HFBSolver",
        "hamiltonian" => dictify(solver.hamiltonian),
        "size" => dictify(solver.size),
        "temperature" => dictify(solver.temperature)
        )
end


function dictify(solution::HFB.HFBSolution)
    OrderedDict(
        "type" => "HFBSolution",
        "rho" => dictify(solution.ρ),
        "t" => dictify(solution.t),
        "Gamma" => dictify(solution.Γ),
        "Delta" => dictify(solution.Δ),
        )
end




objectify(obj::Number) = obj
objectify(obj::AbstractString) = obj

objectify(obj::AbstractArray{T,N}) where {T,N} = begin
    out = [objectify(x) for x in obj]
    DType = Union{}
    for x in out
        DType = typejoin(DType, typeof(x))
    end
    if DType <: Union{}
        DType = Any
    end
    return convert(AbstractArray{DType,N}, out)
end

function objectify(obj::Dict)
    ks = keys(obj)
    return OrderedDict(objectify(k) => objectify(obj[k]) for k in ks)
end

const DICTHANDLER = Dict(
    "Symbol" => ( d -> Symbol(d["value"]) ),
    "Tuple" => ( d -> tuple(objectify(d["value"])...) ),
    "FractCoord" => d -> begin
        whole = objectify(d["whole"])
        fraction = objectify(d["fraction"])
        FractCoord(whole, fraction)
    end,
    "UnitCell" => d -> begin
        latticevectors = objectify(d["latticevectors"])
        orbitals = objectify(d["orbitals"])
        OrbitalType = typejoin([typeof(orb[1]) for orb in orbitals]...)
        OrbitalType = isa(OrbitalType, DataType) ? OrbitalType : Any
        uc = newunitcell(latticevectors; OrbitalType=OrbitalType)
        for orb in orbitals
            addorbital!(uc, orb[1], orb[2])
        end
        return uc
    end,
    "HoppingDiagonal" => d -> begin
        amplitude = objectify(d["amplitude"])
        i = objectify(d["i"])
        Ri = objectify(d["Ri"])
        return Spec.HoppingDiagonal(amplitude, i, Ri)
    end,
    "HoppingOffdiagonal" =>  d -> begin
        amplitude = objectify(d["amplitude"])
        i = objectify(d["i"])
        j = objectify(d["j"])
        Ri = objectify(d["Ri"])
        Rj = objectify(d["Rj"])
        return Spec.HoppingOffdiagonal(amplitude, i, j, Ri, Rj)
    end,
    "InteractionDiagonal" => d -> begin
        amplitude = objectify(d["amplitude"])
        i = objectify(d["i"])
        j = objectify(d["j"])
        Ri = objectify(d["Ri"])
        Rj = objectify(d["Rj"])
        return Spec.InteractionDiagonal(amplitude, i, j, Ri, Rj)
    end,
    "InteractionOffdiagonal" => d-> begin
        amplitude = objectify(d["amplitude"])
        i = objectify(d["i"])
        j = objectify(d["j"])
        k = objectify(d["k"])
        l = objectify(d["l"])
        Ri = objectify(d["Ri"])
        Rj = objectify(d["Rj"])
        Rk = objectify(d["Rk"])
        Rl = objectify(d["Rl"])
        Spec.InteractionOffdiagonal(amplitude,
                                    i, j, k, l,
                                    Ri, Rj, Rk, Rl)
    end,
    "Hamiltonian" => d -> begin
        unitcell = objectify(d["unitcell"])
        hoppings = objectify(d["hoppings"])
        interactions = objectify(d["interactions"])
        hamspec = Spec.Hamiltonian(unitcell)
        for hop in hoppings
            Spec.addhopping!(hamspec, hop)
        end
        for int in interactions
            Spec.addinteraction!(hamspec, int)
        end
        return hamspec
    end,
    "HoppingMeanField" => d -> begin
        amplitude = d["amplitude"]
        target = objectify( d["target"] )
        source = objectify( d["source"] )
        star = d["star"]
        return HFB.HoppingMeanField(amplitude, target, source, star)
    end,

    "PairingMeanField" => d -> begin
        amplitude = d["amplitude"]
        target = objectify( d["target"] )
        source = objectify( d["source"] )
        negate = d["negate"]
        return HFB.PairingMeanField(amplitude, target, source, negate)
    end,

    "HFBHamiltonian" => d -> begin
        unitcell = objectify(d["unitcell"])
        hoppings = Spec.Hopping[objectify(d["hoppings"])...]
        particle_hole_interactions = HFB.HoppingMeanField[objectify(d["particle_hole_interactions"])...]
        particle_particle_interactions = HFB.PairingMeanField[objectify(d["particle_particle_interactions"])...]
        return HFB.HFBHamiltonian(unitcell,
                                  hoppings,
                                  particle_hole_interactions,
                                  particle_particle_interactions)
    end,

    "HFBComputer" => d -> begin
        unitcell = objectify(d["unitcell"])
        hoppings = Spec.Hopping[objectify(d["hoppings"])...]
        temperature = d["temperature"]
        ρ_registry = objectify(d["rho_registry"])
        t_registry = objectify(d["t_registry"])
        Γ_registry = objectify(d["Gamma_registry"])
        Δ_registry = objectify(d["Delta_registry"])
        comp = HFB.HFBComputer(unitcell, hoppings, temperature)
        comp.ρ_registry = ρ_registry
        comp.t_registry = t_registry
        comp.Γ_registry = Γ_registry
        comp.Δ_registry = Δ_registry
    end,

    "HFBSolver" => d -> begin
        hamiltonian = objectify(d["hamiltonian"])
        size = objectify(d["size"])
        temperature = objectify(d["temperature"])
        solver = HFB.HFBSolver(hamiltonian, size, temperature)
        return solver
    end
)


function objectify(dict ::OrderedDict)
    if haskey(DICTHANDLER, dict["type"])
        return DICTHANDLER[dict["type"]](dict)
    else
        return dict
    end
end

end # module Objectify
