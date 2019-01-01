#export hfbfreeenergy
export hfbgrandpotential
export hfbfreeenergy

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

**Should not be used**, since (Γ, Δ) and (ρ, t) are not self-consistent

"""
function hfbfreeenergy_naive(solver::HFBSolver{O},
                             tf::HFBField) where {O}
    fermi = solver.hfbcomputer.fermi
    norb = numorbital(solver.hfbcomputer.unitcell)

    computeH = make_hamiltonian(solver.hfbcomputer, tf)
    computeT = make_hoppingmatrix(solver.hfbcomputer)
    computeΓ = make_Gammamatrix(solver.hfbcomputer, tf)
    computeΔ = make_Deltamatrix(solver.hfbcomputer, tf)

    E_T = zero(Float64)
    E_Γ = zero(Float64)
    E_Δ = zero(Float64)
    S = zero(Float64)
    for k in solver.momentumgrid
        H = computeH(k)
        T = computeT(k)
        Γ = computeΓ(k)
        Δ = computeΔ(k)

        (eigenvalues, eigenvectors) = (eigen(Hermitian(H))...,)
        f = [fermi(e) for e in eigenvalues]
        ψ = reshape(eigenvectors, (2, norb, 2*norb))
        u = ψ[1, :, :]
        v = ψ[2, :, :]

        uf = u * Diagonal(f)
        ρmat = uf * (u')
        tmat = uf * (v')
        for i=1:norb, j=1:norb
            E_T += real( T[i,j] * ρmat[j,i] )
            E_Γ += real( Γ[i,j] * ρmat[j,i] )
            E_Δ += real( Δ[i,j] * conj(tmat[i,j]) )
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


#=
function hopping_energy(hop::HoppingDiagonal{R},
                        ρ ::AbstractMatrix{C},
                        momentum::AbstractVector{K}) where {R <:Real, C <:Number, K <:Real}
                        t, i, Ri = hop.amplitude, hop.i, hop.Ri
    return t * real(ρ[i,i])
end

function hopping_energy(hop::HoppingOffDiagonal{R},
                        ρ ::AbstractMatrix{C},
                        momentum::AbstractVector{K}) where {R <:Real, C <:Number, K <:Real}
    t, i, j, Ri, Rj = hop.amplitude, hop.i, hop.j, hop.Ri, hop.Rj
    return t * real(ρ[i,i])
end
=#

function pairing_energy(int ::InteractionDiagonal{<:Real},
                        t ::AbstractMatrix{<:Number},
                        momentum::AbstractVector{<:Real})
    v, i, j = int.amplitude, int.i, int.j
    return v * real( t[i,j] * conj(t[j,i]) )
end

function pairing_energy(int ::InteractionOffdiagonal{<:Number},
                        t::AbstractMatrix{<:Number},
                        momentum::AbstractVector{<:Real})
    v, i, j, k, l = int.amplitude, int.i, int.j, int.k, int.l
    return 2 * ( v * (tmat[k,l] * conj(t[j,i]) + t[l,k] * conj(t[i,j])) )
end

function hfbfreeenergy(solver::HFBSolver,
                       hfbfield::HFBField;
                       update::Function=simpleupdate)
    fermi = solver.hfbcomputer.fermi
    newhfbfield = let
        newhfbamplitude = next_hfbamplitude(solver, hfbfield)
        newhfbfield = deepcopy(hfbfield)
        update(newhfbfield, make_hfbfield(solver, newhfbamplitude))
        newhfbfield
    end
    return hfbfreeenergy(solver, hfbfield, newhfbfield)
end

function hfbfreeenergy(solver::HFBSolver,
                       tf::HFBField,
                       ntf::HFBField)
    fermi = solver.hfbcomputer.fermi
    norb = numorbital(solver.hamiltonian.unitcell)

    # compute H with old solution, Γ and Δ with new solution
    computeH = make_hamiltonian(solver.hfbcomputer, tf)
    computeT = make_hoppingmatrix(solver.hfbcomputer)
    computeΓ = make_Gammamatrix(solver.hfbcomputer, ntf)
    computeΔ = make_Deltamatrix(solver.hfbcomputer, ntf)

    E_T = zero(Float64)
    E_Γ = zero(Float64)
    E_Δ = zero(Float64)
    S = zero(Float64)
    for k in solver.momentumgrid
        H = computeH(k)
        T = computeT(k)
        Γ = computeΓ(k)
        Δ = computeΔ(k)

        (eigenvalues, eigenvectors) = (eigen(Hermitian(H))...,)
        f = [fermi(e) for e in eigenvalues]
        ψ = reshape(eigenvectors, (2, norb, 2*norb))
        u = ψ[1, :, :]
        v = ψ[2, :, :]

        uf = u * Diagonal(f)
        ρmat = uf * (u')
        tmat = uf * (v')

        for i=1:norb, j=1:norb
            E_T += real( T[i,j] * ρmat[j,i] )
            E_Γ += real( Γ[i,j] * ρmat[j,i] )
            E_Δ += real( Δ[i,j] * conj(tmat[i,j]) )
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
                           solution::HFBAmplitude) where {O}
    E = zero(Float64)
    for (idx_Γ, (diag, i, j, R, srcs)) in enumerate(solver.hfbcomputer.Γ_registry)
        E_Γ = zero(Float64)
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
                           hfbfield::HFBField) where {O}
    fermi = solver.hfbcomputer.fermi
    norb = numorbital(solver.hfbcomputer.unitcell)

    computeH = makehamiltonian(solver.hfbcomputer, hfbfield)
    computeT = makehoppingmatrix(solver.hfbcomputer)
    #computeΓ = makeGammamatrix(solver.hfbcomputer, solution.Γ)
    #computeΔ = makeDeltamatrix(solver.hfbcomputer, solution.Δ)

    E_T = zero(Float64)
    E_V = zero(Float64)
    S = zero(Float64)
    for k in solver.momentumgrid
        H = computeH(k)
        T = computeT(k)
        #Γ = computeΓ(k)
        #Δ = computeΔ(k)

        (eigenvalues, eigenvectors) = (eigen(Hermitian(H))...,)
        f = [fermi(e) for e in eigenvalues]
        ψ = reshape(eigenvectors, (2, norb, 2*norb))
        u = ψ[1, :, :]
        v = ψ[2, :, :]

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
