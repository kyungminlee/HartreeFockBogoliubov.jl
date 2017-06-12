
function fermidirac(temperature ::Real; ttol=eps(Float64), etol=sqrt(eps(Float64)))
  @assert(temperature >= 0)
  if temperature < ttol
    function(e)
      if e < -etol
        return 1.0
      elseif e < etol
        return 0.5
      else
        return 0.0
      end
    end
  else
    beta = 1.0 / temperature
    function(e)
      return 1.0 / (exp(beta * e) + 1.0)
    end
  end
end

#--------------------------- DRAFT

const CollectRow = Tuple{Int64, Int64, Vector{Float64}}
const DeployRow = Tuple{Int64, Int64, Vector{Float64}, Vector{Tuple{Int64, Complex128, Bool}}}

type Solver
  unitcell ::UnitCell
  #size ::Vector{Int64}
  hoppings ::Vector{Embed.Hopping}
  temperature ::Float64
  fermi ::Function
  ρ_registry ::Vector{CollectRow}
  t_registry ::Vector{CollectRow}
  Γ_registry ::Vector{DeployRow}
  Δ_registry ::Vector{DeployRow}
end


function Solver(ham::HFBHamiltonian,
                #size ::AbstractVector{Int64},
                temperature::Real;
                ttol=eps(Float64),
                etol=sqrt(eps(Float64)))
  dim = dimension(ham.unitcell)
  #@assert(length(size) == dim, "size should match dimension")
  #@assert(all((x) -> x>0, size), "size should be all positive")
  @assert(temperature >= 0, "temperature should be non-negative")

  unitcell = ham.unitcell
  hoppings = [Embed.embed(unitcell, hop) for hop in ham.hoppings]
  fermi = fermidirac(temperature)

  function getdistance(i ::Int64, j ::Int64, Rij::AbstractVector{Int64})
    ri, rj = getorbitalcoord(unitcell, i), getorbitalcoord(unitcell, j)
    rj = rj + Rij
    ri, rj = fract2carte(unitcell, ri), fract2carte(unitcell, rj)
    return rj - ri
  end

  collect_reg = Dict()
  deploy_reg = Dict()
  for hopmf in ham.particle_hole_interactions
    let
      (k,l,Rkl) = hopmf.source
      rkl = getdistance(k, l, Rkl)
      collect_reg[k,l,Rkl] = (length(collect_reg)+1, (k, l, rkl))
    end

    let
      v = hopmf.amplitude
      (i,j,Rij) = hopmf.target
      rij = getdistance(i, j, Rij)
      deploy_reg[i,j,Rij] = (length(deploy_reg)+1, (i, j, rij, []))
    end
  end

  for hopmf in ham.particle_hole_interactions
    v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    srcidx = collect_reg[k,l,Rkl][1]
    star = hopmf.targetconj
    push!( deploy_reg[i,j,Rij][2][4], (srcidx, v, star) )
  end

  ρ_registry = sort([(idx, val) for (key, (idx, val)) in collect_reg], by=(x) -> x[1])
  Γ_registry = sort([(idx, val) for (key, (idx, val)) in deploy_reg], by=(x) -> x[1])
  ρ_registry = [val for (idx, val) in ρ_registry]
  Γ_registry = [val for (idx, val) in Γ_registry]

  collect_reg = Dict()
  deploy_reg = Dict()

  for hopmf in ham.particle_particle_interactions
    let
      (k,l,Rkl) = hopmf.source
      rkl = getdistance(k, l, Rkl)
      collect_reg[k,l,Rkl] = (length(collect_reg)+1, (k, l, rkl))
    end

    let
      v = hopmf.amplitude
      (i,j,Rij) = hopmf.target
      rij = getdistance(i, j, Rij)
      deploy_reg[i,j,Rij] = (length(deploy_reg)+1, (i, j, rij, []))
    end
  end

  for hopmf in ham.particle_particle_interactions
    v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    srcidx = collect_reg[(k,l,Rkl)][1]
    neg = hopmf.negate
    push!( deploy_reg[i,j,Rij][2][4], (srcidx, v, neg) )
  end
  t_registry = sort([(idx, val) for (key, (idx, val)) in collect_reg], by=(x) -> x[1])
  Δ_registry = sort([(idx, val) for (key, (idx, val)) in deploy_reg], by=(x) -> x[1])
  t_registry = [val for (idx, val) in t_registry]
  Δ_registry = [val for (idx, val) in Δ_registry]

  return Solver(unitcell,
                #size,
                hoppings,
                temperature, fermi,
                ρ_registry, t_registry,
                Γ_registry, Δ_registry)
end


"""
func : (idx, i, j, r) -> val
"""
function makefields(funcρ ::Function, funct ::Function, solver ::Solver)
  ρs = zeros(Complex128, length(solver.ρ_registry))
  ts = zeros(Complex128, length(solver.t_registry))

  for (idx, (i, j, r)) in enumerate(solver.ρ_registry)
    v = funcρ(idx, i, j, r)
    if i==j && all((x)->x==0, r)
      ρs[idx] = real(v)
    else
      ρs[idx] = v
    end
  end
  for (idx, (i, j, r)) in enumerate(solver.t_registry)
    ts[idx] = funct(idx, i, j, r)
  end
  return (ρs, ts)
end


function makefields(solver ::Solver)
  ρs = zeros(Complex128, length(solver.ρ_registry))
  ts = zeros(Complex128, length(solver.t_registry))
  return (ρs, ts)
end


function computemeanfield(solver ::Solver,
                          ρs ::AbstractVector{Complex128},
                          ts ::AbstractVector{Complex128})
  Γs = zeros(Complex128, length(solver.Γ_registry))
  Δs = zeros(Complex128, length(solver.Δ_registry))
  for (tgtidx, (i, j, r, srcs)) in enumerate(solver.Γ_registry)
    value = 0.0 + 0.00im
    for (srcidx, amplitude, star) in srcs
      value += amplitude * (star ? conj(ρs[srcidx]) : ρs[srcidx])
    end
    Γs[tgtidx] = value
  end

  for (tgtidx, (i, j, r, srcs)) in enumerate(solver.Δ_registry)
    value = 0.0 + 0.00im
    for (srcidx, amplitude, neg) in srcs
      value += amplitude * (neg ? -ts[srcidx] : ts[srcidx])
    end
    Δs[tgtidx] = value
  end
  return (Γs, Δs)
end




function makehamiltonian(solver ::Solver,
                         Γs ::AbstractVector{Complex128},
                         Δs ::AbstractVector{Complex128})
  norb = numorbital(solver.unitcell)
  hk = Generator.generatefast(solver.unitcell, solver.hoppings)
  function(k ::AbstractVector{Float64})
    out = zeros(Complex128, (norb, 2, norb, 2))
    # 1/3. non-interacting kinetic part
    hk( k, view(out, :,1,:,1))
    hk(-k, view(out, :,2,:,2))
    # 2/3. Gamma
    for (idx, Γ) in enumerate(Γs)
      (i, j, r, s) = solver.Γ_registry[idx]
      out[i,1,j,1] += Γ * exp(1im * dot( k, r))
      out[i,2,j,2] += Γ * exp(1im * dot(-k, r))
    end

    out[:,2,:,2] = -transpose(out[:,2,:,2])

    # 3/3. Delta
    for (idx, Δ) in enumerate(Δs)
      (i, j, r, s) = solver.Δ_registry[idx]
      out[i,1,j,2] += Δ * exp(1im * dot( k, r))
      out[j,2,i,1] += conj(Δ) * exp(-1im * dot(k, r))
    end
    return reshape(out, (norb*2, norb*2))
  end
end

function makegreencollectors(solver::Solver)
  fermi = solver.fermi
  norb = numorbital(solver.unitcell)
  ρ_registry = solver.ρ_registry
  t_registry = solver.t_registry

  function(k::AbstractVector{Float64},
           eigenvalues ::AbstractVector{Float64},
           eigenvectors ::AbstractMatrix{Complex128},
           ρout ::AbstractVector{Complex128},
           tout ::AbstractVector{Complex128})
    @assert(length(eigenvalues) == 2*norb)
    @assert(size(eigenvectors) == (2*norb, 2*norb))

    f = [fermi(e) for e in eigenvalues]
    ψ = reshape(eigenvectors, (norb, 2, norb*2))
    u = ψ[:, 1, :]
    v = ψ[:, 2, :]

    function ρfunc(i ::Int64, j ::Int64)
      sum(f .* u[i, :] .* conj(u[j, :]))
    end
    function tfunc(i ::Int64, j ::Int64)
      sum(f .* u[i, :] .* conj(v[j, :]))
    end

    for (idx, (i, j, r)) in enumerate(ρ_registry)
      ρout[idx] += ρfunc(i, j) * exp(-1im * dot(k, r))
    end
    for (idx, (i, j, r)) in enumerate(t_registry)
      tout[idx] += tfunc(i, j) * exp(-1im * dot(k, r))
    end
  end
end


function dumpyaml(solver::Solver)

end


#=
function something(hamiltonian::Hamiltonian, n1 ::Int64, n2 ::Int64)
  uc = hamiltonian.unitcell
  hops = hamiltonian.hoppings
  ph_int = hamiltonian.particle_hole_interactions
  pp_int = hamiltonian.particle_particle_interactions

  ks = generate k points

  Γs = random
  Δs = random

  hk_kinetic = generatehoppingfast(hops)
  hk_hartree = generatehartreefockbogoliubovfast(ph_int, Γs)
  hk_pairing = generatehartreefockbogoliubovfast(pp_int, Δs)

  for k in ks
    hmat = zeros(hamiltonian)
    hk_kinetic(k, hmat)
    hk_hartree(k, hmat)
    hk_pairing(k, hmat)

    @assert is hmat Hermitian?
    hmat = 0.5 * (hmat + conj(transpose(hmat)))
    (eivals, eivecs) = eigs(Hermitian(hmat))
    f = [fermidirac(e) for e in eivals]

    for inter in ph_int
      inter.target
      inter.source
      inter.amplitude
      inter.targetdisplacement
      inter.sourcedisplacement
    end


  end



end
=#



#=

#-----------------------------
type Solver
  unitcell ::UnitCell
  hoppings ::Vector{Embed.Hopping}
  greencollectors
  #ρ_collectors ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
  #t_collectors ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
  #ρ_deployers  ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
  #t_deployers  ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
end



function Solver(model ::Hamiltonian, temperature ::Real)
  lv ::Array{Float64, 2} = model.unitcell.latticevectors
  ph = model.particle_hole_interactions
  pp = model.particle_particle_interactions

  greencollectors = makegreencollectors(model, temperature)

  #ρ_collectors = makecollectors(ph)
  #t_collectors = makecollectors(pp)
  #ρ_deployers  = makedeployers(ph)
  #t_deployers  = makedeployers(pp)

  return Solver(model.unitcell,
                model.hoppings,
                greencollectors,
                ρ_collectors, t_collectors,
                ρ_deployers,  t_deployers)
end



"""
 makegreencollector

Create a function which

 # Returns
"""
function makegreencollectors(hamiltonian::Hamiltonian, temperature ::Float64)

  fermi = fermidirac(temperature)
  norb = numorbital(hamiltonian.unitcell)
  function(eigenvalues ::Vector{Float64}, eigenvectors ::Matrix{Complex128})
    @assert(length(eigenvalues) == 2 * norb)
    @assert(size(eigenvectors) == (2*norb, 2*norb))

    f = [fermi(e) for e in eigenvalues]
    ψ = reshape(eigenvectors, (norb, 2, norb*2))
    u = ψ[:, 1, :]
    v = ψ[:, 2, :]

    function ρfunc(i ::Int64, j ::Int64)
      sum(f .* u[i, :] .* conj(u[j, :]))
    end
    function tfunc(i ::Int64, j ::Int64)
      sum(f .* u[i, :] .* conj(v[j, :]))
    end
    return (ρfunc, tfunc)
  end
end


function makecollectors(hamiltonian::Hamiltonian, greencollectors::Tuple{Any, Any})
  (ρfunc, tfunc) = greencollectors
  unitcell = hamiltonian.unitcell
  ph_int = hamiltonian.particle_hole_interactions
  pp_int = hamiltonian.particle_particle_interactions
  ρcollector = makecollectors(unitcell, ph_int)
  tcollector = makecollectors(unitcell, pp_int)
  #  ρcollector(k, ρfunc)
end


"""
"""
function makecollectors(unitcell::UnitCell, interactions ::Vector{MeanField})
  collectors = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in interactions
    (k, l) = Gamma.source
    rFij::FractCoord = Gamma.targetdisplacement
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCij::CarteCoord = fract2carte(unitcell.lattice, rFij)
    rCkl::CarteCoord = fract2carte(unitcell.lattice, rFkl)

    if haskey(collectors, (k, l, rFkl.w))
      prev_func = collectors[(k, l, rFkl.w)]
      new_func = (momentum::Vector{Float64}, green) -> begin
        return (green(k,l) * exp(-1im * dot(rCkl, momentum))
                + prev_func(momentum, green))
      end
      collectors[(k, l, rFkl.w)] = new_func
    else
      new_func = (momentum::Vector{Float64}, green) -> begin
        return (green(k,l) * exp(-1im * dot(rCkl, momentum)))
      end
      collectors[(k, l, rFkl.w)] = new_func
    end
  end
  return collectors
end


"""
"""
function makedeployers(unitcell::UnitCell, interactions ::Vector{MeanField})
  deployers = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in interactions
    (i, j) = Gamma.target
    (k, l) = Gamma.source
    V = Gamma.amplitude
    rFij::FractCoord = Gamma.targetdisplacement
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCij::CarteCoord = fract2carte(unitcell.lattice, rFij)
    rCkl::CarteCoord = fract2carte(unitcell.lattice, rFkl)

    if haskey(deployers, (i, j, rFij.w))
      prev_func = deployers[(i, j, rFij.w)]
      deployers[(i, j, rFij.w)] = (momentum, values) -> begin
        return (V * values[(k, l, rFkl.w)] * exp(1im * dot(rCij, momentum))
                + prev_func(momentum, values))
      end
    else
      deployers[(i, j, rFij.w)] = (momentum, values) -> begin
        return (V * values[(k, l, rFkl.w)] * exp(1im * dot(rCij, momentum)))
      end
    end
  end
  return deployers
end



function collect(unitcell::UnitCell, ρ_collect, t_collect, momentum, ρ, t)
  ρ_value = Dict()
  t_value = Dict()

  for ((k, l, Rkl), collectfunc) in ρ_collect
    ρ_value[(k, l, Rkl)] = collectfunc(momentum, ρ)
  end

  for ((k, l, Rkl), collectfunc) in t_collect
    t_value[(k, l, Rkl)] = collectfunc(momentum, t)
  end
  return (ρ_value, t_value)
end


function deploy_dense(unitcell ::UnitCell,
                      ρ_deploy,
                      t_deploy,
                      momentum,
                      ρ_value,
                      t_value)
  n = length(unitcell.orbitals)
  Γ = zeros(Complex128, (n, n))
  Δ = zeros(Complex128, (n, n))

  for ((i,j,_), deployfunc) in ρ_deploy
    Γ[i,j] += deployfunc(momentum, ρ_value)
  end
  for ((i,j,_), deployfunc) in t_deploy
    Δ[i,j] += deployfunc(momentum, t_value)
  end
end
=#
