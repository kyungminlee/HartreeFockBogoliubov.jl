
using ..Lattice
using ..Spec
using ..HFB

export nambufy
export freeze

"""
"""
function nambufy(hop::HoppingDiagonal{R}) where {R<:Real}
    hop1 = HoppingDiagonal{R}( hop.amplitude, hop.i*2-1, hop.Ri)
    hop2 = HoppingDiagonal{R}(-hop.amplitude, hop.i*2  , hop.Ri)
    return (hop1, hop2)
end

"""
"""
function nambufy(hop::HoppingOffdiagonal{C}) where {C<:Number}
    hop1 = HoppingOffdiagonal{C}(      hop.amplitude,  hop.i*2-1, hop.j*2-1, hop.Ri, hop.Rj)
    hop2 = HoppingOffdiagonal{C}(-conj(hop.amplitude), hop.i*2  , hop.j*2  , hop.Ri, hop.Rj)
    return (hop1, hop2)
end

"""
"""
function nambufy(unitcell::UnitCell{O}) where {O<:Tuple}
    NewOrbitalType = Tuple{String, O.parameters...}
    nambuunitcell = Lattice.make_unitcell(unitcell.latticevectors; OrbitalType=NewOrbitalType)
    for (orb, fc) in unitcell.orbitals
        addorbital!(nambuunitcell, (orb..., "Particle"), fc)
        addorbital!(nambuunitcell, (orb..., "Hole"), fc)
    end
    return nambuunitcell
end

"""
"""
function nambufy(unitcell::UnitCell{O}) where {O}
    NewOrbitalType = Tuple{String, O}
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
