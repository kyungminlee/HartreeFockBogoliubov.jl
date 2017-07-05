export HFBSolver
export getnextsolution,
       getnextsolutionpython,
       loop,
       simpleupdate

using ProgressMeter


"""
"""
mutable struct HFBSolver{T}
  # Originals
  hamiltonian ::FullHamiltonian{T}
  size ::Vector{Int64}
  temperature ::Float64

  # Derivatives
  momentumgrid ::Array{Vector{Float64}}
  hfbhamiltonian ::HFBHamiltonian{T}
  hfbcomputer ::HFBComputer{T}
  greencollectors ::Function
end


"""
"""
function HFBSolver(hamiltonian::FullHamiltonian{T},
                   size ::AbstractVector{<:Integer},
                   temperature ::Real;
                   tol::Float64=eps(Float64)) where {T}
  dim = dimension(hamiltonian.unitcell)
  @assert(length(size) == dim, "dimension and size do not match: size=$(size)")
  @assert(all((x) -> x > 0, size), "size should be positive: size=$(size)")
  @assert(temperature >= 0, "temperature should be non-negative")

  momentumgrid = Lattice.momentumgrid(hamiltonian.unitcell, size)

  hfbhamiltonian = HFB.HFBHamiltonian(hamiltonian)
  hfbcomputer = HFB.HFBComputer(hfbhamiltonian, temperature)
  greencollectors = HFB.makegreencollectors(hfbcomputer)
  return HFBSolver{T}(hamiltonian, size, temperature,
                      momentumgrid,
                      hfbhamiltonian, hfbcomputer, greencollectors)
end


"""
"""
function getnextsolution(solver ::HFBSolver{T}, sol ::HFBSolution) where {T}
  newsol = newhfbsolution(solver.hfbcomputer)
  ham = makehamiltonian(solver.hfbcomputer, sol.Γ, sol.Δ)
  for momentum in solver.momentumgrid
    hk = ham(momentum)
    (eivals, eivecs) = eig(Hermitian(hk))
    solver.greencollectors(momentum, eivals, eivecs, newsol.ρ, newsol.t)
  end
  newsol.ρ /= length(solver.momentumgrid)
  newsol.t /= length(solver.momentumgrid)
  newsol.Γ[:], newsol.Δ[:] = HFB.computetargetfields(solver.hfbcomputer, newsol.ρ, newsol.t)
  return newsol
end


"""
"""
function getnextsolutionthreaded(solver ::HFBSolver{O},
                                 sol ::HFBSolution) where {O}
  newsol = newhfbsolution(solver.hfbcomputer)
  const ham = makehamiltonian(solver.hfbcomputer, sol.Γ, sol.Δ)

  threadedρ = [zeros(newsol.ρ) for tid in 1:Threads.nthreads()]
  threadedt = [zeros(newsol.t) for tid in 1:Threads.nthreads()]

  num_momentum = length(solver.momentumgrid)
  Threads.@threads for idx_momentum in 1:num_momentum
    momentum = solver.momentumgrid[idx_momentum]
    tid = Threads.threadid()
    hk = ham(momentum)
    (eivals, eivecs) = eig(Hermitian(hk))
    solver.greencollectors(momentum, eivals, eivecs,
                           threadedρ[tid], threadedt[tid])
  end

  for ρ in threadedρ
    newsol.ρ += ρ
  end
  for t in threadedt
    newsol.t += t
  end

  newsol.ρ /= length(solver.momentumgrid)
  newsol.t /= length(solver.momentumgrid)
  newsol.Γ[:], newsol.Δ[:] = HFB.computetargetfields(solver.hfbcomputer, newsol.ρ, newsol.t)
  return newsol
end



function simpleupdate(sol::HFBSolution, newsol ::HFBSolution)
  @assert(length(sol.ρ) == length(newsol.ρ))
  @assert(length(sol.t) == length(newsol.t))
  @assert(length(sol.Γ) == length(newsol.Γ))
  @assert(length(sol.Δ) == length(newsol.Δ))
  sol.ρ[:] = newsol.ρ[:]
  sol.t[:] = newsol.t[:]
  sol.Γ[:] = newsol.Γ[:]
  sol.Δ[:] = newsol.Δ[:]
  sol
end


function loop(solver ::HFBSolver{T},
              sol::HFBSolution,
              run::Integer;
              update::Function=simpleupdate,
              precondition::Function=identity,
              progressbar::Bool=false
              ) where {T}
  sol = copy(sol)
  if progressbar
    @showprogress for i in 1:run
      precondition(sol)
      newsol = getnextsolution(solver, sol)
      update(sol, newsol)
    end
  else
    for i in 1:run
      precondition(sol)
      newsol = getnextsolution(solver, sol)
      update(sol, newsol)
    end
  end
  sol
end


#=
function totalfreeenergy(solver ::HFBSolver{T}, sol::HFBSolution) where T
  error("Not Implemented")
  const computeT = generatehoppingfast(solver.hamiltonian)
  const norb::Int64 = numorbital(solver.unitcell)
  const computeH = makehamiltonian(solver.hfbcomputer, sol.Γ, sol.Δ)

  for k in solver.momentumgrid
    Hmat = computeH(k)
    Tmat = reshape(H, (norb, 2, norb, 2))[:,1,:,1]
    (eigenvalues, eigenvectors) = eig(Hermitian(H))
    ψ = reshape(eigenvectors, (norb, 2, norb*2))
    u = ψ[:, 1, :]
    v = ψ[:, 2, :]
    f = [fermi(e) for e in eigenvalues]

    ρ(i ::Int64, j ::Int64) = sum(f .* u[i, :] .* conj(u[j, :]))
    t(i ::Int64, j ::Int64) = sum(f .* u[i, :] .* conj(v[j, :]))

    for (idx, Γ) in enumerate(sol.Γ)
      (i, j, r, s) = solver.hfbcomputer.Γ_registry[idx]
      Tmat[i, j] += 0.5 * Γ * cis(dot( k, r))
    end

    energyT = 0.0
    for i in 1:norb, j in 1:norb
      energyT += Tmat[i, j] * ρ(j, i)
    end

    energyΔ = 0.0
    for (idx, Δ) in enumerate(sol.Δ)
      (i, j, r, s) = solver.hfbcomputer.Δ_registry[idx]
      if i == j && all(x -> abs(x) < eps(Float64), r)
        energyΔ += 0.5 * Δ * conj( t(i, j) )
      else
        energyΔ += Δ * conj( t(i, j) )
      end
    end
  end
  # Energy = tr ( K + 1/2 Γ⋅ρ  + Δ⋅t + ...)
  # E = Tr [ (T + 1/2 Γ)⋅ρ + 1/2 Δ⋅t†]
  # need to
  # 1. loop over momentum
  # 2. compute T + 1/2 Gamma and Delta,
  # 3. compute eigenvalues and eigenvectors
  # 4. take tensor dot product,
  # 5. sum over result from 4 to compute Energy
  # keep track of f log f + ... for entropy
  # compute E + TS at the end.
  error("NotImplemented")
end
=#
