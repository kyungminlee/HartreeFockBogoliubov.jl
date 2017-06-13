export HFBSolver
export getnextsolution

type HFBSolver{T}
  # Originals
  hamiltonian ::Spec.Hamiltonian{T}
  size ::Vector{Int64}
  temperature ::Float64

  # Derivatives
  momentumgrid ::Array{Vector{Float64}}
  hfbhamiltonian ::HFBHamiltonian{T}
  hfbcomputer ::HFBComputer{T}
  greencollectors ::Function
end


function HFBSolver{T}(hamiltonian::Spec.Hamiltonian{T},
                      size ::Vector{Int64},
                      temperature ::Float64;
                      tol::Float64=eps(Float64))
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




function getnextsolution{T}(solver ::HFBSolver{T}, sol ::HFBSolution)
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
