#export addinteraction!

immutable HoppingMeanField
  amplitude ::Real
  target ::Tuple{Int64, Int64, Vector{Int64}}
  source ::Tuple{Int64, Int64, Vector{Int64}}
  targetconj ::Bool
  sourceconj ::Bool
  function HoppingMeanField(amplitude ::Real,
                            target::Tuple{Int64, Int64, Vector{Int64}},
                            source::Tuple{Int64, Int64, Vector{Int64}})
    (i, j, Rij) = target
    (k, l, Rkl) = source
    targetconj = false
    sourceconj = false
    if i > j
      (i,j) = (j,i)
      Rij = -Rij
      targetconj = !targetconj
    end
    if k > l
      (k,l) = (l,k)
      Rkl = -Rkl
      sourceconj = !sourceconj
    end
    new(amplitude, (i,j,Rij), (k,l,Rkl), targetconj, sourceconj)
  end
end

immutable PairingMeanField
  amplitude ::Real
  target ::Tuple{Int64, Int64, Vector{Int64}}
  source ::Tuple{Int64, Int64, Vector{Int64}}
  negate ::Bool
  function PairingMeanField(amplitude ::Real,
                            target::Tuple{Int64, Int64, Vector{Int64}},
                            source::Tuple{Int64, Int64, Vector{Int64}})
    (i, j, Rij) = target
    (k, l, Rkl) = source
    negate = false
    if i > j
      (i,j) = (j,i)
      Rij = -Rij
      negate = !negate
    end
    if k > l
      (k,l) = (l,k)
      Rkl = -Rkl
      negate = !negate
    end
    new(amplitude, (i,j,Rij), (k,l,Rkl), negate)
  end
end

type HFBHamiltonian
  unitcell ::UnitCell
  hoppings ::Vector{Spec.Hopping}
  particle_hole_interactions ::Vector{HoppingMeanField}
  particle_particle_interactions ::Vector{PairingMeanField}
end


function HFBHamiltonian(full_hamiltonian ::Spec.Hamiltonian)
  unitcell = full_hamiltonian.unitcell
  hoppings = full_hamiltonian.hoppings
  model = HFBHamiltonian(unitcell, hoppings, [], [])
  for interaction in full_hamiltonian.interactions
    addinteraction!(model, interaction)
  end
  return model
end


"""
Add diagonal interaction
"""
function addinteraction!(model ::HFBHamiltonian,
                         specint ::Spec.InteractionDiagonal)
  v = specint.amplitude
  (i, j) = (specint.i,  specint.j )
  (Ri, Rj) = (specint.Ri, specint.Rj)

  Γs = [
    HoppingMeanField(v, (i, i, Ri-Ri), (j, j, Rj-Rj)),
    HoppingMeanField(v, (j, j, Rj-Rj), (i, i, Ri-Ri)),
    HoppingMeanField(v, (j, i, Ri-Rj), (i, j, Rj-Ri)),
  ]

  Δs = [
    PairingMeanField(v, (i, j, Rj-Ri), (i, j, Rj-Ri)),
  ]
  append!(model.particle_hole_interactions, Γs)
  append!(model.particle_particle_interactions, Δs)
end

#=
"""
    addinteraction!

Add offdiagonal interaction
"""
function addinteraction!(model::HFBHamiltonian,
                         embint ::Embed.InteractionOffdiagonal)
  error("Unimplemented")
  V = embint.amplitude
  Vc = conjugate(V)
  (i, j, k, l)  = (embint.i,  embint.j,  embint.k,  embint.l )
  (ri, rj, rk, rl) = (embint.ri, embint.rj, embint.rk, embint.rl)

  # Hermitian -> completely expanded out
  Γs = [
    MeanField((i,k), (j,l),       V, ri - rk, rj - rl),
    MeanField((j,k), (i,l),      -V, rj - rk, ri - rl),
    MeanField((i,l), (j,k),      -V, ri - rl, rj - rk),
    MeanField((j,l), (i,k),       V, rj - rl, ri - rk),
    # Hermitian conjugate
    MeanField((k,i), (l,j),      Vc, rk - ri, rl - rj),
    MeanField((k,j), (l,i),     -Vc, rk - rj, rl - ri),
    MeanField((l,i), (k,j),     -Vc, rl - ri, rk - rj),
    MeanField((l,j), (k,i),      Vc, rl - rj, rk - ri),
  ]
  Δs = [
    MeanField((i,j), (k,l),   0.5*V, ri - rj, rk - rl),
    MeanField((i,j), (l,k),  -0.5*V, ri - rj, rl - rk),
    MeanField((j,i), (k,l),  -0.5*V, rj - ri, rk - rl),
    MeanField((j,i), (l,k),   0.5*V, rj - ri, rl - rk),
    # Hermitian conjugate
    MeanField((k,l), (i,j),  0.5*Vc, rk - rl, ri - rj),
    MeanField((k,l), (j,i), -0.5*Vc, rk - rl, rj - ri),
    MeanField((l,k), (i,j), -0.5*Vc, rl - rk, ri - rj),
    MeanField((l,k), (j,i),  0.5*Vc, rl - rk, rj - ri),
  ]
  append!(model.particle_hole_interactions, Γs)
  append!(model.particle_particle_interactions, Δs)
end
=#
