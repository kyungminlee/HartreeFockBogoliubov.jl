#export addinteraction!


immutable MeanField
  target ::Tuple{Int64, Int64}
  source ::Tuple{Int64, Int64}
  amplitude ::Real
  targetdisplacement ::Vector{Float64}
  sourcedisplacement ::Vector{Float64}
end


type Hamiltonian
  unitcell ::UnitCell
  hoppings ::Vector{Embed.Hopping}
  particle_hole_interactions ::Vector{MeanField}
  particle_particle_interactions ::Vector{MeanField}
end


function Hamiltonian(full_hamiltonian ::Embed.Hamiltonian)
  unitcell = full_hamiltonian.unitcell
  hoppings = full_hamiltonian.hoppings
  model = Hamiltonian(unitcell, hoppings, [], [])
  for interaction in full_hamiltonian.interactions
    addinteraction!(model, interaction)
  end
  return model
end


function addinteraction!(model ::Hamiltonian,
                         embint ::Embed.InteractionDiagonal)
  V = embint.amplitude
  (i, j)  = (embint.i,  embint.j )
  ri = fract2carte(model.unitcell, embint.ri)
  rj = fract2carte(model.unitcell, embint.rj)

  Γ = [
    MeanField((i,i), (j,j),      V, ri - ri, rj - rj),
    MeanField((j,i), (i,j),     -V, rj - ri, ri - rj),
    MeanField((i,j), (j,i),     -V, ri - rj, rj - ri),
    MeanField((j,j), (i,i),      V, rj - rj, ri - ri),
  ]
  Δ = [
    MeanField((i,j), (i,j),  0.5*V, ri - rj, ri - rj),
    MeanField((i,j), (j,i), -0.5*V, ri - rj, rj - ri),
    MeanField((j,i), (i,j), -0.5*V, rj - ri, ri - rj),
    MeanField((j,i), (j,i),  0.5*V, rj - ri, rj - ri),
  ]
  append!(model.particle_hole_interactions, Γ)
  append!(model.particle_particle_interactions, Δ)
end


function addinteraction!(model::Hamiltonian,
                         embint ::Embed.InteractionOffdiagonal)
  error("Unimplemented")
  V = embint.amplitude
  Vc = conjugate(V)
  (i, j, k, l)  = (embint.i,  embint.j,  embint.k,  embint.l )
  (ri,rj,rk,rl) = (embint.ri, embint.rj, embint.rk, embint.rl)

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
