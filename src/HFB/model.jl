#export addinteraction!


immutable DiagonalHoppingSubscript
  i ::Int64
end

# i < j
"""
required i < j
rij = rj - ri
"""
immutable OffdiagonalHoppingSubscript
  i ::Int64
  j ::Int64
  Rij ::Vector{Int64}
  star ::Bool
end


# i < j
immutable OffdiagonalPairingSubscript
  i ::Int64
  j ::Int64
  Rij ::Vector{Int64}
  negative ::Bool
end

const HoppingSubscript = Union{DiagonalHoppingSubscript,
                               OffdiagonalHoppingSubscript}

const PairingSubscript = Union{OffdiagonalPairingSubscript}


function hoppingsubscript(i ::Int64, j ::Int64, Rij ::Vector{Int64}, star::Bool=false)
  if i == j && all((x) -> (x==0), Rij)
    return DiagonalHoppingSubscript(i)
  elseif i > j
    return OffdiagonalHoppingSubscript(j, i, -Rij, !star)
  else
    return OffdiagonalHoppingSubscript(i, j, Rij, star)
  end
end


function pairingsubscript(i ::Int64, j ::Int64, Rij ::Vector{Int64}, negative::Bool=false)
  @assert(!(i == j && isapprox(rij, zeros(Float64, length(rij)))),
          "pairing vanishes!")
  if i == j && all((x) -> (x==0), Rij)
    @assert false
  elseif i > j
    return OffdiagonalPairingSubscript(j, i, -Rij, !negative)
  else
    return OffdiagonalPairingSubscript(i, j, Rij, negative)
  end
end


immutable HoppingMeanField
  target ::HoppingSubscript
  source ::HoppingSubscript
  amplitude ::Real
end


immutable PairingMeanField
  target ::PairingSubscript
  source ::PairingSubscript
  amplitude ::Real
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
  V = specint.amplitude
  (i, j) = (specint.i,  specint.j )
  (Ri, Rj) = (specint.Ri, specint.Rj)

  Γs = let
    HS = hoppingsubscript
    R0 = zeros(Ri)
    [
      HoppingMeanField(HS(i, i, Ri-Ri), HS(j, j, Rj-Rj), V),
      HoppingMeanField(HS(j, j, Rj-Rj), HS(i, i, Ri-Ri), V),
      HoppingMeanField(HS(j, i, Ri-Rj), HS(i, j, Rj-Ri), V),
    ]
  end

  Δs = let
    PS = pairingsubscript
    [
      PairingMeanField(PS(i, j, Rj-Ri), PS(i, j, Rj-Ri), V),
    ]
  end
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
