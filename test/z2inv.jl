#=
  Computing Z2 invariant
  ======================

  1. Construct (Kane-Mele-Hubbard) Hamiltonian
  2. Construct a hopping-hamiltonian-generator H1
  3. Squarify the Lattice and make generator H2
  4. Make time-reversal-accounted Index Grid and Momentum Grid 
  5. Define time-reversal operation unitary matrix U_T
     5-1. Confirm unitarity
     5-2. Confirm U_T^2 = -1
     5-3. Confirm h(k) = [ U_T ⋅ h( (-k) mod G ) ⋅ U_T^(-1) ]^* for all k
     5-4. Confirm h(k) = [ U_T ⋅ h(k) ⋅ U_T^(-1) ]^*
  6. Contruct eigenvector grid? (select numbers of them)?
     for TRI, POSINT, POS, (and NEG from POS)
  7. Compute F for POSINT & POS/TRI/NEG with ky = 0
  8. Compute A for POS, NEG, TRI on both ky = 0 and ky = π
=#
using DataStructures
import PyPlot


using HartreeFockBogoliubov
import HartreeFockBogoliubov: Spec, Generator, HFB, Topology
using HartreeFockBogoliubov: HFB

function makekanemelehamiltonian(μ ::Float64,
                                 t ::Float64,
                                 mAB ::Float64,
                                 λIsing ::Float64,
                                 U ::Float64,
                                 V ::Float64)
  a0 = [ 0.0, 0.0]
  a1 = [ 0.0, 1.0]
  a2 = [-sqrt(3.0) * 0.5,-0.5]
  a3 = [ sqrt(3.0) * 0.5,-0.5]
  b1 = a3 - a2
  b2 = a1 - a3
  b3 = a2 - a1

  #c1 = b1 - b2
  #c2 = b2 - b3

  unitcell = newunitcell([b1 b2], OrbitalType=Tuple{Symbol, Symbol})
  for spin in [:UP, :DN]
    addorbital!(unitcell, (:A, spin), carte2fract(unitcell, a0))
    addorbital!(unitcell, (:B, spin), carte2fract(unitcell, a1))
  end

  hamspec = Spec.FullHamiltonian(unitcell)

  orbloc = Dict(
    :A => a0,
    :B => a1,
  )

  for sp in [:UP, :DN]
    let orb = :A
      r = orbloc[orb]
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ + mAB, (orb, sp), r))
    end
    let orb = :B
      r = orbloc[orb]
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ - mAB, (orb, sp), r))
    end
  end

  for sp in [:UP, :DN]
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, sp), (:B, sp), a0, a1))
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, sp), (:B, sp), a0, a2))
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, (:A, sp), (:B, sp), a0, a3))
  end

  if abs(λIsing) > eps(Float64)
    for sp in [:UP, :DN], orb in [:A, :B]
      λ = 1im * λIsing * (sp == :UP ? 1 : -1) * (orb == :A ? 1 : -1)
      r = orbloc[orb]
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, (orb, sp), (orb, sp), r, r + b1))
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, (orb, sp), (orb, sp), r, r + b2))
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -λ, (orb, sp), (orb, sp), r, r + b3))
    end
  end

  if abs(U) > eps(Float64)
    for orb in [:A, :B]
      r = orbloc[orb, idx]
      Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, U, (orb, :UP), (orb, :DN), r, r))
    end
  end

  if abs(V) > eps(Float64)
    for sp1 in [:UP, :DN], sp2 in [:UP, :DN]
      Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, sp1), (:B, sp2), a0, a1))
      Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, sp1), (:B, sp2), a0, a2))
      Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, (:A, sp1), (:B, sp2), a0, a3))
    end
  end
  return hamspec
end


function spinhalftimereversal(uc::UnitCell{O}, spinindex) where {O<:Tuple}
  @assert(
    all(let 
      sp = orbname[spinindex]
      sp == :UP || sp == :DN
    end for (orbname, fc) in uc.orbitals) 
  )

  rows = Int[]
  cols = Int[]
  vals = Float64[]
  for (idx, (orbname, fc)) in enumerate(uc.orbitals)
    sp = orbname[spinindex]
    (newsp, sgn) = (sp == :UP) ? (:DN, 1) : (:UP, -1)
    sgn = (sp == :UP) ? 1 : -1
    neworbname = ((i != spinindex ? v : newsp for (i, v) in enumerate(orbname))...)
    jdx = getorbitalindex(uc, neworbname)
    push!(rows, idx)
    push!(cols, jdx)
    push!(vals, sgn)
  end
  return sparse(rows, cols, vals)
end


function main1()
  Eg = 0.4
  for mAB in linspace(0, Eg, 11)
    λIsing = (Eg - mAB) / (3.0*sqrt(3.0))
    kmh1 = makekanemelehamiltonian(0.0, 1.0, mAB, λIsing, 0.0, 0.0)
    #@show numorbital(kmh1.unitcell)
    #@show kmh1.unitcell
    timereversalmatrix = spinhalftimereversal(kmh1.unitcell, 2)
    #@show timereversalmatrix

    z2index = Topology.z2invariant(
                kmh1.unitcell,
                kmh1.hoppings,
                timereversalmatrix,
                8,
                8,
                1:1)
    @printf("%f\t%d\n", mAB, z2index)
  end
end

#=
function main2()
  kmh1 = makekanemelehamiltonian(0.0, 1.0, 0.0, 0.2, 0.0, 0.0)
  kmh2 = Topology.squarify(kmh1)
  norb = numorbital(kmh1.unitcell)
  ham1 = Generator.generatehoppingfast(kmh1)
  ham2 = Generator.generatehoppingfast(kmh2)

  kg1 = momentumgrid(kmh1.unitcell, [8, 8])
  kg2 = momentumgrid(kmh1.unitcell, [8, 8])

  #for k in kg
  #  out = zeros(Complex128, (norb, norb))
  #  hk = ham(k, out)
  #end
  #@show out

  @show kmh2.unitcell.orbitals
  timereversalmatrix = let
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for (idx, (orbname, fc)) in enumerate(kmh2.unitcell.orbitals)
      sp = orbname[end]
      (newsp, sgn) = (sp == :UP) ? (:DN, 1) : (:UP, -1)
      sgn = (sp == :UP) ? 1 : -1
      neworbname = (orbname[1:end-1]..., newsp)
      jdx = getorbitalindex(kmh2.unitcell, neworbname)
      push!(rows, idx) # convention?
      push!(cols, jdx)
      push!(vals, sgn)
      println("$idx -> $jdx")
    end
    sparse(rows, cols, vals)
  end
  @show timereversalmatrix
  @show full(timereversalmatrix)

  g = Topology.triindexgrid([4,3])
  kg2 = momentumgrid(kmh2.unitcell, [8,6])
  for (pc1, (pt, pc2)) in g
    k = kg2[(pc1+1)...]
    if pt == :TRI
      hk = zeros(Complex128, (norb, norb))
      #hkprime = zeros(Complex128, (norb, norb))
      ham2(k, hk)
      #ham2(-k, hkprime)
      hkprime = timereversalmatrix * hk * timereversalmatrix'
      @show k
      #@show hk
      #@show hkprime
      @show maximum(abs.(hk - hkprime))
    elseif pt == :POS
    elseif pt == :NEG
    elseif pt == :POSINT
    elseif pt == :NEGINT
    end
  end

  return
  g = Topology.triindexgrid([4,3])
  fig = PyPlot.figure()
  for (pc1, (pt, pc2)) in g
    if pt == :TRI
      PyPlot.plot(pc1[1], pc1[2], "ro", alpha=0.8)
    elseif pt == :POS
      PyPlot.plot(pc1[1], pc1[2], "go", alpha=0.8)
    elseif pt == :NEG
      PyPlot.plot(pc1[1], pc1[2], "bo", alpha=0.8)
    elseif pt == :POSINT
      PyPlot.plot(pc1[1], pc1[2], "go", alpha=0.2)
    elseif pt == :NEGINT
      PyPlot.plot(pc1[1], pc1[2], "bo", alpha=0.2)
    end
  end
  fig[:savefig]("foo.png")
end
=#

main1()