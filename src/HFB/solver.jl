
function fermidirac(temperature ::Real; ttol=eps(Float64), etol=sqrt(eps(Float64)))
  @assert(temperature >= 0)
  if temperature < ttol
    function(e)
      if e < -tol
        return 1.0
      elseif e < tol
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


type Solver
  unitcell ::UnitCell
  size ::Vector{Int64}
  hoppings ::Vector{Embed.Hopping}
  fermi
  ρ_registry
  t_registry
  Γ_registry
  Δ_registry
end


function Solver(ham::HFBHamiltonian,
                size ::Vector{Int64},
                temperature::Real;
                ttol=eps(Float64),
                etol=sqrt(eps(Float64)))
  dim = dimension(ham.unitcell)
  @assert(length(size) == dim, "size should match dimension")
  @assert(all((x) -> x>0, size), "size should be all positive")
  @assert(temperature >= 0, "temperature should be non-negative")

  unitcell = ham.unitcell
  hoppings = ham.hoppings
  fermi = fermidirac(temperature)


  for interaction in ham.particle_hole_interactions
    tgt = interaction.target
    src = interaction.source
    V = interaction.amplitude
    collect_reg[src] = (length(collect_reg) + 1, rCkl, false)
  end


  function register(interactions)
    collect_reg = Dict()
    deploy_reg = Dict()
    for interaction in interactions
      tgt = interaction.target
      src = interaction.source
      V = interaction.amplitude
      rFij = tgt.rij
      rFkl = tgt.rkl
      rCij = fract2carte(unitcell, rFij)
      rCkl = fract2carte(unitcell, rFkl)

      if !( (k, l, rFkl.whole) in keys(collect_reg) )
        collect_reg[k, l, rFkl.whole] = (length(collect_reg) + 1, rCkl, false)
      end
      if !( (i, j, rFij.whole) in keys(deploy_reg) )
        deploy_reg[i, j, rFij.whole] = (length(deploy_reg) + 1, rCij, [])
      end
    end

    for interaction in interactions
      (i,j) = interaction.target
      (k,l) = interaction.source
      V = interaction.amplitude
      rFij = interaction.targetdisplacement
      rFkl = interaction.sourcedisplacement
      rCij = fract2carte(unitcell, rFij)
      rCkl = fract2carte(unitcell, rFkl)

      collect_idx = collect_reg[k, l, rFkl.whole]
      push!( deploy_reg[i, j, rFij.whole][3], (V, collect_idx) )
    end

    collect_reg = [(idx, (i, j, r)) for ((i, j, R), (idx, r)) in collect_reg]
    deploy_reg = [(idx, (i, j, r, source)) for ((i, j, R), (idx, r, source)) in deploy_reg]
    collect_registry = [val for (idx, val) in sort(collect_reg)]
    deploy_registry = [val for (idx, val) in sort(deploy_reg)]
    return (collect_registry, deploy_registry)
  end

  (ρ_registry, Γ_registry) = register(ham.particle_hole_interactions)
  (t_registry, Δ_registry) = register(ham.particle_particle_interactions)

  #=
  println("ρ_registry")
  println("i\tj\tr")
  for (i, j, r) in ρ_registry
    println("$(i)\t$(j)\t$(r)")
  end
  println()
  println("Γ_registry")
  println("i\tj\tr")
  for (i, j, r, sources) in Γ_registry
    println("$(i)\t$(j)\t$(r):")
    for s in sources
      println("  $(s)")
    end
  end

  println("t_registry")
  println("i\tj\tr")
  for (i, j, r) in t_registry
    println("$(i)\t$(j)\t$(r)")
  end
  println()
  println("Δ_registry")
  println("i\tj\tr")
  for (i, j, r, sources) in Δ_registry
    println("$(i)\t$(j)\t$(r):")
    for s in sources
      println("  $(s)")
    end
  end
  =#
  return Solver(unitcell, size,
                hoppings, fermi,
                ρ_registry, t_registry,
                Γ_registry, Δ_registry)
end


"""
func : (idx, i, j, r) -> val
"""
function makefields(funcΓ, funcΔ, solver ::Solver)
  Γ = zeros(Complex128, length(solver.Γ_registry))
  Δ = zeros(Complex128, length(solver.Δ_registry))

  for (idx, (i, j, r, s)) in enumerate(solver.Γ_registry)
    Γ[idx] = funcΓ(idx, i, j, r)
  end
  for (idx, (i, j, r, s)) in enumerate(solver.Δ_registry)
    Δ[idx] = funcΔ(idx, i, j, r)
  end
  return (Γ, Δ)
end


function makehamiltonian(solver ::Solver,
                         Γs ::Vector{Complex128},
                         Δs ::Vector{Complex128})
  norb = numorbital(solver.unitcell)
  hk = Generator.generatefast(solver.unitcell, solver.hoppings)
  function(k ::Vector{Float64})
    out = zeros(Complex128, (norb, 2, norb, 2))
    # 1/3. non-interacting kinetic part
    hk(k,  view(out, :,1,:,1))
    hk(-k, view(out, :,2,:,2))
    # 2/3. Gamma
    for (idx, Γ) in enumerate(Γs)
      (i, j, r, s) = solver.Γ_registry[idx]
      out[i,1,j,1] += Γ * exp(1im * dot(k, r))
      out[i,2,j,2] += Γ * exp(1im * dot(-k, r))
    end

    out[:,2,:,2] = -transpose(out[:,2,:,2])

    # 3/3. Delta
    for (idx, Δ) in enumerate(Δs)
      (i, j, r, s) = solver.Δ_registry[idx]
      out[i,1,j,2] += Δ * exp(1im * dot(k, r))
      out[j,2,i,1] += conj(Δ) * exp(-1im * dot(k, r))
    end
    return reshape(out, (2*norb, 2*norb))
  end
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
