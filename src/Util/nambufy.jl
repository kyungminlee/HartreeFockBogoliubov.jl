
using ..Lattice
using ..Spec
using ..HFB

export nambufy
export freeze

"""
    nambufy

Makes two copies of the given hopping element, one in the particle basis, and another in the hole basis.
The Nambu space index is the least significant subindex of the new basis.
In other words, given a hopping from `i` to `j` with amplitude `t` (i and j starting from 1), the function `nambufy` returns
two hoppings, one from `2*i-1` to `2*j-1` with amplitude `t`, and another from `2*i` to `2*j` with amplitude `-conj(t)`.
"""
function nambufy(hop::HoppingDiagonal{R}) where {R<:Real}
    hop1 = HoppingDiagonal{R}( hop.amplitude, hop.i*2-1, hop.Ri)
    hop2 = HoppingDiagonal{R}(-hop.amplitude, hop.i*2  , hop.Ri)
    return (hop1, hop2)
end

function nambufy(hop::HoppingOffdiagonal{C}) where {C<:Number}
    hop1 = HoppingOffdiagonal{C}(      hop.amplitude,  hop.i*2-1, hop.j*2-1, hop.Ri, hop.Rj)
    hop2 = HoppingOffdiagonal{C}(-conj(hop.amplitude), hop.i*2  , hop.j*2  , hop.Ri, hop.Rj)
    return (hop1, hop2)
end

"""
    nambufy

Generate a "nambufied" version of a given unitcell, i.e. double the number of orbitals by considering the "particle"
orbitals and the "hole" orbitals. The new orbitals acquire a new name which follows the row-first convention
of Julia. For `UnitCell{O}`, if `O` is a tuple of type `Tuple{A, B, ...}`, the new orbital type is
`Tuple{A, B, ..., String}`. If `O` is not `Tuple`, then the new orbital type is `Tuple{O, String}`.

For example, given a `UnitCell` with orbitals `["UP", "DN"]`, the new `UnitCell` has orbitals
`[("UP", "Particle"), ("UP", "Hole"), ("DN", "Particle"), ("DN", "Hole")]`
If, on the other hand, the original `UnitCell` has orbitals
```julia
[
    ("A", "UP"),
    ("A", "DN"),
    ("B", "UP"),
    ("B", "DN"),
]
```
then the new `UnitCell` has orbitals
```julia
[
    ("A", "UP", "Particle"),
    ("A", "UP", "Hole"),
    ("A", "DN", "Particle"),
    ("A", "DN", "Hole"),
    ("B", "UP", "Particle"),
    ("B", "UP", "Hole"),
    ("B", "DN", "Particle"),
    ("B", "DN", "Hole"),
]
```
"""
function nambufy(unitcell::UnitCell{O}) where {O<:Tuple}
    NewOrbitalType = Tuple{O.parameters..., String}
    nambuunitcell = Lattice.make_unitcell(unitcell.latticevectors; OrbitalType=NewOrbitalType)
    for (orb, fc) in unitcell.orbitals
        addorbital!(nambuunitcell, (orb..., "Particle"), fc)
        addorbital!(nambuunitcell, (orb..., "Hole"), fc)
    end
    return nambuunitcell
end

function nambufy(unitcell::UnitCell{O}) where {O}
    NewOrbitalType = Tuple{O, String}
    nambuunitcell = Lattice.make_unitcell(unitcell.latticevectors; OrbitalType=NewOrbitalType)
    for (orb, fc) in unitcell.orbitals
        addorbital!(nambuunitcell, (orb, "Particle"), fc)
        addorbital!(nambuunitcell, (orb, "Hole"), fc)
    end
    return nambuunitcell
end


"""
    freeze

    Create a hopping hamiltonian out of HFB Hamiltonian and Solution

    Order:
    electrons......, HOLE......
"""
function freeze(computer ::HFBComputer{O},
                hfbfield ::HFBField) where {O}
    norb = numorbital(computer.unitcell)
    unitcell = computer.unitcell
    nambuunitcell = nambufy(unitcell)
    hoppings_diagonal = HoppingDiagonal[]
    hoppings_offdiagonal = HoppingOffiagonal[]

    for hop in computer.hoppings_diagonal
        (hop1, hop2) = nambufy(hop)
        push!(hoppings_diagonal, hop1)
        push!(hoppings_diagonal, hop2)
    end
    for hop in computer.hoppings_offdiagonal
        (hop1, hop2) = nambufy(hop)
        push!(hoppings_offdiagonal, hop1)
        push!(hoppings_offdiagonal, hop2)
    end

    # 2/3. Gamma
    for (idx, Γ) in enumerate(tf.Γ)
        (isdiag, i, j, r, _) = computer.Γ_registry[idx]
        if isdiag
            @assert(i==j && all(x -> isapprox(x, 0.0), r))
            @assert(isapprox(imag(Γ), 0.0))

            iname, icoord = getorbital(unitcell, i)
            Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))

            push!(hoppings_diagonal, HoppingDiagonal{Float64}( real(Γ), i*2-1, Ri))
            push!(hoppings_diagonal, HoppingDiagonal{Float64}(-real(Γ), i*2  , Ri))
        else
            @assert(!(i==j && all(x -> isapprox(x, 0.0), r)))

            iname, icoord = getorbital(unitcell, i)
            jname, jcoord = getorbital(unitcell, j)
            Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
            Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)

            push!(hoppings_offdiagonal, HoppingOffdiagonal{ComplexF64}( Γ      , i*2-1, j*2-1, Ri, Rj))
            push!(hoppings_offdiagonal, HoppingOffdiagonal{ComplexF64}(-conj(Γ), i*2  , j*2  , Ri, Rj))
        end
    end

    # 3/3. Delta
    for (idx, Δ) in enumerate(tf.Δ)
        #info("Adding $Δ to the frozen hoppings")
        (isdiag, i, j, r, _) = computer.Δ_registry[idx]

        @assert( !isdiag )
        @assert( !(i==j && all(x -> isapprox(x, 0.0), r)) )

        iname, icoord = getorbital(unitcell, i)
        jname, jcoord = getorbital(unitcell, j)
        Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
        Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)

        #@show ( Δ, i, j + norb, Ri, Rj)
        #@show (-Δ, j, i + norb, Rj, Ri)
        push!(hoppings_offdiagonal, HoppingOffdiagonal{ComplexF64}( Δ , i*2-1, j*2, Ri, Rj))
        push!(hoppings_offdiagonal, HoppingOffdiagonal{ComplexF64}(-Δ , j*2-1, i*2, Rj, Ri))
    end
    return (nambuunitcell, hoppings_diagonal, hoppings_offdiagonal)
end
