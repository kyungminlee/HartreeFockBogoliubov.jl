type HFBSolver
  # Originals
  hamiltonian ::Spec.Hamiltonian
  size ::Vector{Int64}
  temperature ::Float64

  # Derivatives
  momentumgrid ::Vector{Vector{Float64}}
  hfbhamiltonian ::HFBHamiltonian
  hfbcomputer ::HFBComputer
  greencollectors ::Tuple{Function, Function}

  function HFBSolver(hamiltonian::Spec.Hamiltonian,
                     size ::Vector{Int64},
                     temperature ::Float64;
                     tol::Float64=eps(Float64))
    dim = dimension(hamiltonian.unitcell)
    @assert(length(size) == dim, "dimension and size do not match: size=$(size)")
    @assert(all((x) -> x > 0), "size should be positive: size=$(size)")
    @assert(temperature >= 0, "temperature should be non-negative")

    momentumgrid = Lattice.momentumgrid(hamiltonian.unitcell, size)

    hfbhamiltonian = HFB.HFBHamiltonian(hamiltonian)
    hfbcomputer = HFB.HFBComputer(hfbhamiltonian, temperature)
    greencollectors = HFB.makegreencollectors(hfbcomputer)
    return HFBSolver(hamiltonian, size, temperature,
                     momentumgrid,
                     hfbhamiltonian, hfbcomputer, greencollectors)
  end
end

type HFBSolution
  ρ ::Vector{Complex128}
  t ::Vector{Complex128}
  Γ ::Vector{Complex128}
  Δ ::Vector{Complex128}
end

function newhfbsolution(solver::HFBSolver)
  ρ, t = HFB.makesourcefields(solver.computer)
  Γ, Δ = HFB.computetargetfields(computer, ρ, t)
  solution = HFBSolution(ρ, t, Γ, Δ)
end


function randomize!(sol ::HFBSolution, amplitude=1.0)
  sol.ρ[:] = (rand(Float64, length(sol.ρ)) .* 2 .- 1) * amplitude
  sol.t[:] = (rand(Complex128, length(sol.t)) .* 2 .- 1) * amplitude
  sol.Γ[:], sol.Δ[:] = HFB.computetargetfields(computer, sol.ρ, sol.t)
end


function next(solver ::HFBSolver, sol ::HFBSolution)
  newsol = newhfbsolution(solver)
  ham = makehamiltonian(solver.hfbcomputer, sol.Γ, sol.Δ)
  for momentum in solver.momentumgrid
    hk = ham(k)
    (eivals, eivecs) = eig(Hermitian(hk))
    solver.greencollectors(k, eivals, eivecs, newsol.ρ, newsol.t)
  end
  newsol.ρ /= length(solver.momentumgrid)
  newsol.t /= length(solver.momentumgrid)
  newsol.Γ[:], newsol.Δ[:] = HFB.computetargetfields(computer, newsol.ρ, newsol.t)
  return newsol
end
