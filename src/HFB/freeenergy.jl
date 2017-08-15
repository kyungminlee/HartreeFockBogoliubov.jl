#export hfbfreeenergy
export hfbgrantpotential


raw"""

Compute Ω = E - T S
where
```math
E = \mathrm{Tr}\left[ (T + \frac{1}{2} \Gamma) \rho + \frac{1}{2} \Delta t^{\dagger} \right]
```
```math
S = - \sum_{i} f_i \log f_i
```
Here i runs through all nambu indices.

"""
function hfbfreeenergy(solver::HFBSolver{O},
                       solution::HFBSolution) where {O}
  fermi = solver.hfbcomputer.fermi
  norb = numorbital(solver.hfbcomputer.unitcell)

  computeH = makehamiltonian(solver.hfbcomputer, solution.Γ, solution.Δ)
  computeT = makehoppingmatrix(solver.hfbcomputer)
  computeΓ = makeGammamatrix(solver.hfbcomputer, solution.Γ)
  computeΔ = makeDeltamatrix(solver.hfbcomputer, solution.Δ)

  E_T = 0.0
  E_Γ = 0.0
  E_Δ = 0.0
  S = 0.0
  for k in solver.momentumgrid
    H = computeH(k)
    T = computeT(k)
    Γ = computeΓ(k)
    Δ = computeΔ(k)

    (eigenvalues, eigenvectors) = eig(Hermitian(H))
    f = [fermi(e) for e in eigenvalues]
    ψ = reshape(eigenvectors, (norb, 2, norb*2))
    u = ψ[:, 1, :]
    v = ψ[:, 2, :]

    ρ(i::Int, j::Int) = sum(f .* u[i, :] .* conj(u[j, :]))
    t(i::Int, j::Int) = sum(f .* u[i, :] .* conj(v[j, :]))

    for i=1:norb, j=1:norb
      E_T += real( T[i,j] * ρ(j,i) )
      E_Γ += real( Γ[i,j] * ρ(j,i) )
      E_Δ += real( Δ[i,j] * conj(t(i,j)) )
    end
    # f already contains both positive and negative energy states
    fp = [max(x, sqrt(eps(Float64))) for x in f]
    S -= sum( fp .* log.(fp) )
  end
  E = E_T + 0.5 * E_Γ + 0.5 * E_Δ
  F = E - solver.hfbcomputer.temperature * S
  Ω = F / length(solver.momentumgrid)
  return (E,S,Ω)
end

export hfbfreeenergy


function hoppingenergy(hopping::Spec.HoppingDiagonal{R},
                       ρ::Function) where {R<:Real}
  v = hopping.amplitude
  i = hopping.i
  return v * real( ρ(i,i) )
end

function hoppingenergy(hopping::Spec.HoppingOffdiagonal{C},
                       ρ::Function) where {C<:Number}

end




# Faulty
function interactionenergy(interaction::Spec.InteractionDiagonal{R},
                           ρ::Function,
                           t::Function)where {R<:Real}
  v = interaction.amplitude
  i = interaction.i
  j = interaction.j
  val = v * ( real(ρ(i,i)) * real(ρ(j,j)) - abs2(ρ(i,j)) + abs2(t(i,j)))
end


# OK?
function interactionenergy(solver::HFBSolver{O},
                           solution::HFBSolution) where {O}
  E = 0.0
  for (idx_Γ, (diag, i, j, R, srcs)) in enumerate(solver.hfbcomputer.Γ_registry)
    E_Γ = 0.0
    idx_ρ0 = 0
    for (idx_ρ2, (diag2, i2, j2, R2)) in enumerate(solver.hfbcomputer.ρ_registry)
      if i == i2 && j == j2 && R == R2
        idx_ρ0 = idx_ρ2
      end
    end
    ρ0 = solution.ρ[idx_ρ0]
    for (idx_ρ, amplitude, star) in srcs
      E_Γ += amplitude * solution.ρ[idx_ρ]
    end
    E_Γ = real( ρ0 * E_Γ )
    E += E_Γ
  end
  return E
end


function hfbgrandpotential(solver::HFBSolver{O},
                           solution::HFBSolution) where {O}
  fermi = solver.hfbcomputer.fermi
  norb = numorbital(solver.hfbcomputer.unitcell)

  computeH = makehamiltonian(solver.hfbcomputer, solution.Γ, solution.Δ)
  computeT = makehoppingmatrix(solver.hfbcomputer)
  #computeΓ = makeGammamatrix(solver.hfbcomputer, solution.Γ)
  #computeΔ = makeDeltamatrix(solver.hfbcomputer, solution.Δ)

  E_T = 0.0
  E_V = 0.0
  S = 0.0
  for k in solver.momentumgrid
    H = computeH(k)
    T = computeT(k)
    #Γ = computeΓ(k)
    #Δ = computeΔ(k)

    (eigenvalues, eigenvectors) = eig(Hermitian(H))
    f = [fermi(e) for e in eigenvalues]
    ψ = reshape(eigenvectors, (norb, 2, norb*2))
    u = ψ[:, 1, :]
    v = ψ[:, 2, :]

    ρ(i::Int, j::Int) = sum(f .* u[i, :] .* conj(u[j, :]))
    t(i::Int, j::Int) = sum(f .* u[i, :] .* conj(v[j, :]))
    for i=1:norb, j=1:norb
      E_T += real( T[i,j] * ρ(j,i) )
    end
    for interaction in solver.hamiltonian.interactions
      E_V += interactionenergy(interaction, ρ, t)
    end
    # f already contains both positive and negative energy states
    fp = [max(x, sqrt(eps(Float64))) for x in f]
    S -= sum( fp .* log.(fp) )
  end
  E = E_T + E_V
  F = E - solver.hfbcomputer.temperature * S
  Ω = F / length(solver.momentumgrid)
  return (E,S,Ω)
end
