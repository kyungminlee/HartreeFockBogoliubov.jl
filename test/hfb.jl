
#=
@testset "single-site" begin

  latticevectors = [1.0 0.0; 0.0 1.0]
  unitcell = newunitcell(latticevectors, OrbitalType=Symbol)

  addorbital!(unitcell, :UP, FractCoord([0, 0], [0.0, 0.0]))
  addorbital!(unitcell, :DN, FractCoord([0, 0], [0.0, 0.0]))

  spec_hamiltonian = Spec.Hamiltonian(unitcell)
  Spec.addhopping!(spec_hamiltonian,
                   Spec.hoppingbycarte(unitcell, 0.1, :UP, [0.0, 0.0]))
  Spec.addhopping!(spec_hamiltonian,
                   Spec.hoppingbycarte(unitcell,-0.1, :DN, [0.0, 0.0]))

  Spec.addinteraction!(spec_hamiltonian,
      Spec.interactionbycarte(unitcell, -4.0, :UP, :DN, [0.0, 0.0], [0.0, 0.0])
  )

  hfb_hamiltonian = HFB.HFBHamiltonian(spec_hamiltonian)
  @show hfb_hamiltonian

  computer = HFB.HFBComputer(hfb_hamiltonian, [1,1], 0)

  @show computer

  for i in 1:length(computer.ρ_registry)
    println("ρs[$i] = 1.0")
    (ρs, ts) = HFB.makesourcefields(computer)
    ρs[i] = 1.0
    @show ρs
    @show ts

    (Γs, Δs) = HFB.computetargetfields(computer, ρs, ts)
    @show Γs
    @show Δs
  end
  println()
  for i in 1:length(computer.t_registry)
    println("ts[$i] = 1.0")
    (ρs, ts) = HFB.makesourcefields(computer)
    ts[i] = 1.0
    @show ρs
    @show ts

    (Γs, Δs) = HFB.computetargetfields(computer, ρs, ts)
    @show Γs
    @show Δs
  end

  greencollectors = HFB.makegreencollectors(computer)
  ρs, ts = HFB.makesourcefields((x...) -> rand(Float64), (x...) -> rand(Float64), computer)
  (Γs, Δs) = HFB.computetargetfields(computer, ρs, ts)

  ρn, tn = HFB.makesourcefields(computer)

  for run in 1:100
    @show run
    ham = HFB.makehamiltonian(computer, Γs, Δs)
    ρn[:] = 0
    tn[:] = 0
    k = [0.0, 0.0]
    hk = ham(k)
    #hk = 0.5 * (hk + conj(transpose(hk)))
    @show hk
    (eivals, eivecs) = eig(Hermitian(hk))
    greencollectors(k, eivals, eivecs, ρn, tn)

    (Γn, Δn) = HFB.computetargetfields(computer, ρn, tn)

    @show ρn
    @show tn
    @show Γn
    @show Δn
    (Γs, Δs) = (Γn * 0.5 + Γs * 0.5, Δn * 0.5 + Δs * 0.5)
  end





  #=
  computer = HFB.HFBComputer(hfb_hamiltonian, [1,1], 0)



  let
    hk = Generator.generatefast(computer.unitcell, computer.hoppings)
    out = zeros(Complex128, (2,2))
    hk([0.0, 0.0], out)
  end

  println()
  println()
  println()
  println()



  rng = MersenneTwister(1)

  (Γs, Δs) = HartreeFockBogoliubovModel.makesourcefields(
                (x...) -> rand(rng, Float64),
                (x...) -> rand(rng, Complex128),
                 computer)

  for idx in 1:length(Γs)
    println("setting Γ[$(idx)] = 100.0")
    Γs[:] = 0.0
    Δs[:] = 0.0
    Γs[idx] = 100.0
    hmfk = HartreeFockBogoliubovModel.makehamiltonian(computer, Γs, Δs)
    hmat = hmfk([0.0, 0.0])

    for i in 1:4
      for j in 1:4
        print("$(hmat[i,j])\t")
      end
      println()
    end
  end
  # TODO: HERE: MEAN FIELD and Hamiltonian, and etc.
  # 1. Naming.
  # 2. Decide on FractCoord vs CarteCoord
  # 3. Implement HFBComputer
  # 4. Test
  =#

end
=#

@testset "squarehubbard" begin
  latticevectors = [1.0 0.0 ; 0.0 1.0]
  unitcell = newunitcell(latticevectors, OrbitalType=Symbol)
  addorbital!(unitcell, :UP, FractCoord([0,0], [0.0, 0.0]))
  addorbital!(unitcell, :DN, FractCoord([0,0], [0.0, 0.0]))

  for U in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0]
    hamspec = Spec.Hamiltonian(unitcell)
    μ = 0.0
    #U = 5.0
    for sp in [:UP, :DN]
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ - 0.5*U, sp, [0.0, 0.0]))
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -1.0, sp, sp, [0.0, 0.0], [ 1.0, 0.0]))
      Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -1.0, sp, sp, [0.0, 0.0], [ 0.0, 1.0]))
    end

    Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, U, :UP, :DN, [0.0, 0.0], [0.0, 0.0]))
    hamhfb = HFB.HFBHamiltonian(hamspec)
    computer = HFB.HFBComputer(hamhfb, [16, 16], 0.1)
    greencollectors = HFB.makegreencollectors(computer)

    ρs, ts = HFB.makesourcefields((x...) -> rand(Float64), (x...) -> rand(Float64), computer)
    ρn, tn = HFB.makesourcefields(computer)

    #@show computer.ρ_registry
    let
      ρ0 = 0.5 * (ρs[1] + ρs[2])
      ρz = 0.5 * (ρs[1] - ρs[2])
      ρp = ρs[3]

      newρz = sqrt(abs(ρz)^2 + abs(ρp)^2)
      ρs[:] = [ρ0 + newρz, ρ0 - newρz, 0.0+0.0im]
      ts[:] = 0.0
    end

    Γs, Δs = HFB.computetargetfields(computer, ρs, ts)

    kxs = linspace(0, 2*π, 16+1)[1:end-1]
    kys = linspace(0, 2*π, 16+1)[1:end-1]

    for run in 1:100
      ham = HFB.makehamiltonian(computer, Γs, Δs)
      ρn[:] = 0
      tn[:] = 0
      for kx in kxs, ky in kys
        k = [kx, ky]
        hk = ham(k)
        (eivals, eivecs) = eig(Hermitian(hk))
        greencollectors(k, eivals, eivecs, ρn, tn)
      end
      ρn /= (16*16)
      tn /= (16*16)

      Γn, Δn = HFB.computetargetfields(computer, ρn, tn)
      Γdiff = Γn - Γs
      Γs = 0.5 * (Γs + Γn)
      #if maximum(abs.(Γdiff)) < 1E-10
      #   break
      #end
      #@show run
      #@show ρn
      #@show Γn
    end
    println("U = $(U), ρ = $(ρn)")
  end
end
