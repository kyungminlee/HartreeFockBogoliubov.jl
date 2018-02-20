using Base.Test

using DataStructures

using HartreeFockBogoliubov
import HartreeFockBogoliubov: Spec, Generator, HFB
using HartreeFockBogoliubov: HFB

"""
maketriplethamiltonian

# Arguments
* `μ ::Float64`
* `t ::Float64`
* `mAB ::Float64`
* `λIsing ::Float64`
* `U ::Float64`
* `V ::Float64`
"""
function maketripletmdhamiltonian(μ ::Float64,
                                  t ::Float64,
                                  mAB ::Float64,
                                  λIsing ::Float64,
                                  U ::Float64, V ::Float64)
    a0 = [ 0.0, 0.0]
    a1 = [ 0.0, 1.0]
    a2 = [-sqrt(3.0) * 0.5,-0.5]
    a3 = [ sqrt(3.0) * 0.5,-0.5]
    b1 = a3 - a2
    b2 = a1 - a3
    b3 = a2 - a1

    c1 = b1 - b2
    c2 = b2 - b3

    unitcell = newunitcell([c1 c2], OrbitalType=Tuple{Symbol, Int64, Symbol})
    for spin in [:UP, :DN]
        addorbital!(unitcell, (:A, 1, spin), carte2fract(unitcell, a0))
        addorbital!(unitcell, (:B, 1, spin), carte2fract(unitcell, a1))
        addorbital!(unitcell, (:A, 2, spin), carte2fract(unitcell, b1))
        addorbital!(unitcell, (:B, 2, spin), carte2fract(unitcell, b1+a1))
        addorbital!(unitcell, (:A, 3, spin), carte2fract(unitcell, b1+b2))
        addorbital!(unitcell, (:B, 3, spin), carte2fract(unitcell, b1+b2+a1))
    end

    hamspec = Spec.FullHamiltonian(unitcell)

    orbloc = Dict(
            (:A, 1) => a0,
            (:A, 2) => b1,
            (:A, 3) => b1+b2,
            (:B, 1) => a1,
            (:B, 2) => b1+a1,
            (:B, 3) => b1+b2+a1,
    )

    for sp in [:UP, :DN], idx in [1,2,3]
        let orb = :A
            r = orbloc[orb, idx]
            Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ + mAB, (orb, idx, sp), r))
        end
        let orb = :B
            r = orbloc[orb, idx]
            Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ - mAB, (orb, idx, sp), r))
        end
    end

    for sp in [:UP, :DN]
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 1, sp), (:B, 1, sp), a0, a1))
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 1, sp), (:B, 2, sp), a0, a2))
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 1, sp), (:B, 3, sp), a0, a3))

        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 2, sp), (:B, 2, sp), b1, b1+a1))
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 2, sp), (:B, 3, sp), b1, b1+a2))
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 2, sp), (:B, 1, sp), b1, b1+a3))

        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 3, sp), (:B, 3, sp), b1+b2, b1+b2+a1))
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 3, sp), (:B, 1, sp), b1+b2, b1+b2+a2))
        Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, 3, sp), (:B, 2, sp), b1+b2, b1+b2+a3))
    end

    if abs(λIsing) > eps(Float64)
        for sp in [:UP, :DN], idx in [1,2,3], orb in [:A, :B]
            λ = 1im * λIsing * (sp == :UP ? 1 : -1) * (orb == :A ? 1 : -1)
            jdx = mod(idx, 3) + 1
            r = orbloc[orb, idx]
            Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, (orb, idx, sp), (orb, jdx, sp), r, r + b1))
            Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, (orb, idx, sp), (orb, jdx, sp), r, r + b2))
            Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, (orb, idx, sp), (orb, jdx, sp), r, r + b3))
        end
    end

    if abs(U) > eps(Float64)
        for idx in [1,2,3], orb in [:A, :B]
            r = orbloc[orb, idx]
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, U, (orb, idx, :UP), (orb, idx, :DN), r, r))
        end
    end

    if abs(V) > eps(Float64)
        for sp1 in [:UP, :DN], sp2 in [:UP, :DN]
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 1, sp1), (:B, 1, sp2), a0, a1))
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 1, sp1), (:B, 2, sp2), a0, a2))
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 1, sp1), (:B, 3, sp2), a0, a3))

            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 2, sp1), (:B, 2, sp2), b1, b1 + a1))
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 2, sp1), (:B, 3, sp2), b1, b1 + a2))
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 2, sp1), (:B, 1, sp2), b1, b1 + a3))

            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 3, sp1), (:B, 3, sp2), b1 + b2, b1+b2+a1))
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 3, sp1), (:B, 1, sp2), b1 + b2, b1+b2+a2))
            Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, 3, sp1), (:B, 2, sp2), b1 + b2, b1+b2+a3))
        end
    end
    return hamspec
end


@testset "freezetest" begin
    μ = 0.0
    t = 1.0
    mAB = 0.0
    λIsing = 0.2
    U = -3.0
    V = 0.0
    hamspec = maketripletmdhamiltonian(μ, t, mAB, λIsing, U, V)
    hfbsolver = HFBSolver(hamspec, [12, 12], 0.0)

    hfbsolution = newhfbsolution(hfbsolver.hfbcomputer)
    randomize!(hfbsolver.hfbcomputer, hfbsolution)

    hamiltonian1 = makehamiltonian(hfbsolver, hfbsolution.Γ, hfbsolution.Δ)

    (nambuunitcell, nambuhoppings) = HFB.freeze(hfbsolver, hfbsolution.Γ, hfbsolution.Δ)

    hamiltonian2 = let
        foo = Generator.generatefast(nambuunitcell, nambuhoppings)
        norb = numorbital(nambuunitcell)
        function ret(k ::AbstractVector{Float64}) ::Matrix{Complex128}
            out = zeros(Complex128, (norb, norb))
            foo(k, out)
            return out
        end
    end


    for i in 1:100
        k = rand(Float64, 2)*20-10
        @show maximum(abs.(hamiltonian1(k) - hamiltonian2(k)))
    end

end
