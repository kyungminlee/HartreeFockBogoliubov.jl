using DataStructures
import PyPlot


using HartreeFockBogoliubov
using HartreeFockBogoliubov.Spec
using HartreeFockBogoliubov.Generator
using HartreeFockBogoliubov.Topology
using HartreeFockBogoliubov.HFB

function makehaldanemodel(μ ::Float64,
                          t ::Float64,
                          mAB ::Float64,
                          λIsing ::Float64,
                          V ::Float64)
  a0 = [ 0.0, 0.0]
  a1 = [ 0.0, 1.0]
  a2 = [-sqrt(3.0) * 0.5,-0.5]
  a3 = [ sqrt(3.0) * 0.5,-0.5]
  b1 = a3 - a2
  b2 = a1 - a3
  b3 = a2 - a1

  unitcell = make_unitcell([b1 b2], OrbitalType=Symbol)
  addorbital!(unitcell, :A, carte2fract(unitcell, a0))
  addorbital!(unitcell, :B, carte2fract(unitcell, a1))

  hamspec = Spec.FullHamiltonian(unitcell)

  orbloc = Dict(
    :A => a0,
    :B => a1,
  )

  let orb = :A
    r = orbloc[orb]
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ + mAB, orb, r))
  end
  let orb = :B
    r = orbloc[orb]
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ - mAB, orb, r))
  end

  Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, :A, :B, a0, a1))
  Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, :A, :B, a0, a2))
  Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, :A, :B, a0, a3))

  if abs(λIsing) > eps(Float64)
    for orb in [:A, :B]
      λ = 1im * λIsing * (orb == :A ? 1 : -1)
      r = orbloc[orb]
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, orb, orb, r, r + b1))
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, orb, orb, r, r + b2))
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, orb, orb, r, r + b3))
    end
  end

  if abs(V) > eps(Float64)
    Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, :A, :B, a0, a1))
    Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, :A, :B, a0, a2))
    Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, :A, :B, a0, a3))
  end
  return hamspec
end


function main()
  haldanemodel = makehaldanemodel(0.0, 1.0, 0.0, 0.5, 1.0)
  #=
  ham = Generator.generatehoppingfast(haldanemodel)
  norb = numorbital(haldanemodel.unitcell)
  kg = momentumgrid(haldanemodel.unitcell, [16, 16])

  for (i, k) in enumerate(kg)
    out = zeros(ComplexF64, (norb, norb))
    ham(k, out)

  end
  =#
  @show Topology.chernnumber(haldanemodel.unitcell, haldanemodel.hoppings, 16, 16, 1:1)

  hfb_haldanemodel = HFBHamiltonian(haldanemodel)
  hfb_computer = HFBComputer(hfb_haldanemodel, 0.0)
  hfb_hfbamplitude = make_hfbamplitude(hfb_computer)

  srand(1)
  randomize!(hfb_computer, hfb_hfbamplitude)
  hfb_hfbfield = make_hfbfield(hfb_computer, hfb_hfbamplitude)

  @show hfb_hfbfield

  hkgen1 = make_hamiltonian(hfb_computer, hfb_hfbfield)
  (nambuunitcell, ham2_diagonal, ham2_offdiagonal) = freeze(hfb_computer, hfb_hfbfield)

  hkgen2 = let
    norb = numorbital(haldanemodel.unitcell)*2
    generator = Generator.hopping_inplace(nambuunitcell, ham2_diagonal, ham2_offdiagonal)
    function(k::AbstractVector{<:Real})
      out = zeros(ComplexF64, (norb, norb))
      generator(k, out)
      return out
    end
  end

  kgrid = momentumgrid(haldanemodel.unitcell, [3,4])
  for k in kgrid
    hk1 = hkgen1(k)
    hk2 = hkgen2(k)
    @show k
    @show hk1
    @show hk2
    @show maximum(abs.(hk1 - hk2))
    println()
  end
end


main()
