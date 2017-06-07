@testset "single-site" begin
  latticevectors = [1.0 0.0; 0.0 1.0]
  unitcell = newunitcell(latticevectors, OrbitalType=Symbol)

  addorbital!(unitcell, :UP, FractCoord([0, 0], [0.0, 0.0]))
  addorbital!(unitcell, :DN, FractCoord([0, 0], [0.0, 0.0]))

  spec_hamiltonian = Spec.Hamiltonian(unitcell)
  Spec.addhopping!(spec_hamiltonian,
                   Spec.hoppingbycarte(unitcell,  1.0, :UP, [0.0, 0.0]))
  Spec.addhopping!(spec_hamiltonian,
                   Spec.hoppingbycarte(unitcell, -1.0, :DN, [0.0, 0.0]))

  Spec.addinteraction!(spec_hamiltonian,
      Spec.interactionbycarte(unitcell, 5.0, :UP, :DN, [0.0, 0.0], [0.0, 0.0])
  )

  embed_hamiltonian = Embed.Hamiltonian(spec_hamiltonian)
  hfb_hamiltonian = HartreeFockBogoliubovModel.Hamiltonian(embed_hamiltonian)
  @show hfb_hamiltonian

  @assert(false, "READ TODO")
  # TODO: HERE: MEAN FIELD and Hamiltonian, and etc.
  # 1. Naming.
  # 2. Decide on FractCoord vs CarteCoord
  # 3. Implement Solver
  # 4. Test


end
