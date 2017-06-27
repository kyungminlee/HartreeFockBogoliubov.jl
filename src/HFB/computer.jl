export fermidirac
export HFBComputer
export HFBSolution
export HFBHint
export makesourcefields,
       computetargetfields,
       makehamiltonian,
       makegreencollectors,
       newhfbhint,
       newhfbsolution,
       randomize!,
       fixhfbsolution


"""
    fermidirac

Return the Fermi-Dirac distribution function for the given temperature.
For energy whose absolute value is less than `etol`, return 0.5.
"""
function fermidirac(temperature ::T;
                    ttol::T=eps(T),
                    etol::T=sqrt(eps(T))) where {T <:AbstractFloat}
  @assert(temperature >= 0)
  @assert(ttol >= 0)
  @assert(etol >= 0)
  if temperature < ttol
    function(e ::T)
      if e <= -etol
        return T(1.0)
      elseif e < etol
        return T(0.5)
      else
        return T(0.0)
      end
    end
  else
    beta = 1.0 / temperature
    function(e ::T)
      return T(1.0 / (exp(beta * e) + 1.0))
    end
  end
end


"""
    fermidirac

Return the Fermi-Dirac distribution function for the given temperature.
For energy whose absolute value is less than `etol`, return 0.5.
"""
function fermidirac(temperature ::T;
                    ttol::Float64=eps(Float64),
                    etol::Float64=sqrt(eps(Float64))) where {T <:Integer}
  @assert(temperature >= 0, "temperature should be non-negative")
  @assert(ttol >= 0, "ttol should be non-negative")
  @assert(etol >= 0, "etol should be non-negative")
  if temperature == 0
    function(e ::Float64)
      if e <= -etol
        return 1.0
      elseif e < etol
        return 0.5
      else
        return 0.0
      end
    end
  else
    beta = 1.0 / Float64(temperature)
    function(e ::Float64)
      return 1.0 / (exp(beta * e) + 1.0)
    end
  end
end


"""
`CollectRow` is holds info on how to compute ρ or t.
Its elements are:
1. row orbital
2. col orbital
3. displacement r(col) - r(row)
"""
const CollectRow = Tuple{Int64, Int64, Vector{Float64}}


"""
`DeployRow` is holds info on how to compute Γ or Δ.
Its elements are:
1. row orbital
2. col orbital
3. displacement r(col) - r(row)
4. list of sources, each of which is a tuple of
  1. index of ρ or t from which to compute this Γ or Δ.
  2. amplitude (coefficient to multiply to ρ or t)
  3. boolean indicating whether
     (1) conjugation is needed (for ρ/Γ) or
     (2) minus sign is needed (for t/Δ).
"""
const DeployRow = Tuple{Int64, Int64, Vector{Float64}, Vector{Tuple{Int64, Complex128, Bool}}}


"""
`HFBConmputer` is a type holding the ρ, t and Γ, Δ of a Hartree-Fock-Bogoliubov
Hamiltonian.
"""
mutable struct HFBComputer{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{Spec.SpecHopping}
  temperature ::Float64

  fermi ::Function
  ρ_registry ::Vector{CollectRow}
  t_registry ::Vector{CollectRow}
  Γ_registry ::Vector{DeployRow}
  Δ_registry ::Vector{DeployRow}
end

function HFBComputer(unitcell::UnitCell{O},
                     hoppings::AbstractVector{Spec.SpecHopping},
                     temperature::Real) where {O}
 fermi = fermidirac(temperature)
 return HFBComputer{O}(unitcell, hoppings, temperature, fermi, [], [], [], [])
end

function HFBComputer(ham::HFBHamiltonian{T},
                     temperature::Real;
                     ttol=eps(Float64),
                     etol=sqrt(eps(Float64))) where {T}
  @assert(ttol >= 0.0, "ttol should be non-negative")
  @assert(etol >= 0.0, "etol should be non-negative")
  @assert(temperature >= 0, "temperature should be non-negative")
  dim = dimension(ham.unitcell)

  unitcell = ham.unitcell
  hoppings = ham.hoppings
  fermi = fermidirac(temperature)

  function getdistance(i ::Int64, j ::Int64, Rij::AbstractVector{Int64})
    ri, rj = getorbitalcoord(unitcell, i), getorbitalcoord(unitcell, j)
    rj = rj + Rij
    ri, rj = fract2carte(unitcell, ri), fract2carte(unitcell, rj)
    return rj - ri
  end

  collect_reg = Dict()
  deploy_reg = Dict()

  # always collect density
  let
    R = zeros(Int64, dimension(ham.unitcell))
    r = zeros(Float64, dimension(ham.unitcell))
    for (i, orb) in enumerate(ham.unitcell.orbitals)
      collect_reg[i, i, R] = (length(collect_reg)+1, (i, i, r))
    end
  end

  let
    R = zeros(Int64, dimension(ham.unitcell))
    r = zeros(Float64, dimension(ham.unitcell))
    for (i, orb) in enumerate(ham.unitcell.orbitals)
      deploy_reg[i, i, R] = (length(deploy_reg)+1, (i, i, r, []))
    end
  end

  for hopmf in ham.particle_hole_interactions
    let
      (k,l,Rkl) = hopmf.source
      rkl = getdistance(k, l, Rkl)
      if !haskey(collect_reg, (k,l, Rkl))
        collect_reg[k,l,Rkl] = (length(collect_reg)+1, (k, l, rkl))
      end
    end

    let
      v = hopmf.amplitude
      (i,j,Rij) = hopmf.target
      rij = getdistance(i, j, Rij)
      if !haskey(deploy_reg, (i,j, Rij))
        deploy_reg[i,j,Rij] = (length(deploy_reg)+1, (i, j, rij, []))
      end
    end
  end

  for hopmf in ham.particle_hole_interactions
    v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    srcidx = collect_reg[k,l,Rkl][1]
    star = hopmf.star
    if abs(v) > eps(Float64)
      push!( deploy_reg[i,j,Rij][2][4], (srcidx, v, star) )
    end
  end

  let indices = [idx for (key, (idx, val)) in collect_reg]
    @assert(length(unique(indices)) == length(indices))
  end
  let indices = [idx for (key, (idx, val)) in deploy_reg]
    @assert(length(unique(indices)) == length(indices))
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
      if !haskey(collect_reg, (k,l,Rkl))
        collect_reg[k,l,Rkl] = (length(collect_reg)+1, (k, l, rkl))
      end
    end

    let
      v = hopmf.amplitude
      (i,j,Rij) = hopmf.target
      rij = getdistance(i, j, Rij)
      if !haskey(deploy_reg, (i,j,Rij))
        deploy_reg[i,j,Rij] = (length(deploy_reg)+1, (i, j, rij, []))
      end
    end
  end

  for hopmf in ham.particle_particle_interactions
    v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    srcidx = collect_reg[(k,l,Rkl)][1]
    neg = hopmf.negate
    if abs(v) > eps(Float64)
      push!( deploy_reg[i,j,Rij][2][4], (srcidx, v, neg) )
    end
  end

  let indices = [idx for (key, (idx, val)) in collect_reg]
    @assert(length(unique(indices)) == length(indices))
  end
  let indices = [idx for (key, (idx, val)) in deploy_reg]
    @assert(length(unique(indices)) == length(indices))
  end

  t_registry = sort([(idx, val) for (key, (idx, val)) in collect_reg], by=(x) -> x[1])
  Δ_registry = sort([(idx, val) for (key, (idx, val)) in deploy_reg], by=(x) -> x[1])
  t_registry = [val for (idx, val) in t_registry]
  Δ_registry = [val for (idx, val) in Δ_registry]

  return HFBComputer{T}(unitcell,
                        hoppings,
                        temperature, fermi,
                        ρ_registry, t_registry,
                        Γ_registry, Δ_registry)
end


mutable struct HFBSolution
  ρ ::Vector{Complex128}
  t ::Vector{Complex128}
  Γ ::Vector{Complex128}
  Δ ::Vector{Complex128}
end

mutable struct HFBHint{T}
  ρ ::Dict{Tuple{T, T, Vector{Int64}}, Complex128}
  t ::Dict{Tuple{T, T, Vector{Int64}}, Complex128}
end

import Base: copy

copy(x::HFBSolution) = HFBSolution(copy(x.ρ), copy(x.t), copy(x.Γ), copy(x.Δ))
copy{T}(x::HFBHint{T}) = HFBHint{T}(copy(x.ρ), copy(x.t))

import Base: -

function -(sol1::HFBSolution, sol2::HFBSolution)
  HFBSolution(sol1.ρ - sol2.ρ, sol1.t - sol2.t,
              sol1.Γ - sol2.Γ, sol1.Δ - sol2.Δ)
end


"""
func : (idx, i, j, r) -> val
"""
function makesourcefields(funcρ ::Function,
                          funct ::Function,
                          computer ::HFBComputer{O}) where {O}
  ρs = zeros(Complex128, length(computer.ρ_registry))
  ts = zeros(Complex128, length(computer.t_registry))

  for (idx, (i, j, r)) in enumerate(computer.ρ_registry)
    v = funcρ(idx, i, j, r)
    if i==j && all((x)->x==0, r)
      ρs[idx] = real(v)
    else
      ρs[idx] = v
    end
  end
  for (idx, (i, j, r)) in enumerate(computer.t_registry)
    ts[idx] = funct(idx, i, j, r)
  end
  return (ρs, ts)
end


function makesourcefields(computer ::HFBComputer{O}) where {O}
  ρs = zeros(Complex128, length(computer.ρ_registry))
  ts = zeros(Complex128, length(computer.t_registry))
  return (ρs, ts)
end


function computetargetfields(computer ::HFBComputer{O},
                             ρs ::AbstractVector{Complex128},
                             ts ::AbstractVector{Complex128}) where {O}
  Γs = zeros(Complex128, length(computer.Γ_registry))
  Δs = zeros(Complex128, length(computer.Δ_registry))
  for (tgtidx, (i, j, r, srcs)) in enumerate(computer.Γ_registry)
    value = 0.0 + 0.00im
    for (srcidx, amplitude, star) in srcs
      value += amplitude * (star ? conj(ρs[srcidx]) : ρs[srcidx])
    end
    Γs[tgtidx] = value
  end

  for (tgtidx, (i, j, r, srcs)) in enumerate(computer.Δ_registry)
    value = 0.0 + 0.00im
    for (srcidx, amplitude, neg) in srcs
      value += amplitude * (neg ? -ts[srcidx] : ts[srcidx])
    end
    Δs[tgtidx] = value
  end
  return (Γs, Δs)
end


"""
"""
function makehamiltonian(computer ::HFBComputer,
                         Γs ::AbstractVector{Complex128},
                         Δs ::AbstractVector{Complex128})
  norb = numorbital(computer.unitcell)
  hk = Generator.generatefast(computer.unitcell, computer.hoppings)
  function(k ::AbstractVector{Float64})
    out = zeros(Complex128, (norb, 2, norb, 2))
    # 1/3. non-interacting kinetic part
    hk( k, view(out, :,1,:,1))
    hk(-k, view(out, :,2,:,2))
    # 2/3. Gamma
    for (idx, Γ) in enumerate(Γs)
      (i, j, r, s) = computer.Γ_registry[idx]
      out[i,1,j,1] += Γ * exp(1im * dot( k, r))
      out[i,2,j,2] += Γ * exp(1im * dot(-k, r))
    end

    out[:,2,:,2] = -transpose(out[:,2,:,2])

    # 3/3. Delta
    for (idx, Δ) in enumerate(Δs)
      (i, j, r, s) = computer.Δ_registry[idx]
      out[i,1,j,2] += Δ * exp(1im * dot( k, r))
      out[j,2,i,1] += conj(Δ) * exp(-1im * dot(k, r))
    end
    return reshape(out, (norb*2, norb*2))
  end
end


"""
    makegreencollectors

Returns a function which has the following signature
```
collector(k, eigenvalues, eigenvectors, ρout, tout)
```
"""
function makegreencollectors(computer::HFBComputer{O}) where {O}
  fermi = computer.fermi
  norb = numorbital(computer.unitcell)
  ρ_registry = computer.ρ_registry
  t_registry = computer.t_registry

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


function newhfbsolution(computer::HFBComputer{O}) where {O}
  ρ, t = HFB.makesourcefields(computer)
  Γ, Δ = HFB.computetargetfields(computer, ρ, t)
  return HFBSolution(ρ, t, Γ, Δ)
end


"""
    Check if hint contains ρ
"""
function newhfbsolution(computer::HFBComputer{O}, hint::HFBHint{O}) where {O}
  ρ, t = HFB.makesourcefields(computer)
  unitcell = computer.unitcell

  for (ρidx, (i, j, r)) in enumerate(computer.ρ_registry)
    iname, icoord = getorbital(unitcell, i)
    jname, jcoord = getorbital(unitcell, j)
    Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
    Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)
    R = Rj - Ri
    if haskey(hint.ρ, (iname, jname, R))
      ρ[ρidx] = hint.ρ[iname, jname, R]
    elseif haskey(hint.ρ, (jname, iname, -R))
      ρ[ρidx] = conj(hint.ρ[jname, iname, -R])
    end
  end

  for (tidx, (i, j, r)) in enumerate(computer.t_registry)
    iname, icoord = getorbital(unitcell, i)
    jname, jcoord = getorbital(unitcell, j)
    Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
    Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)
    R = Rj - Ri
    if haskey(hint.t, (iname, jname, R))
      t[tidx] = hint.t[iname, jname, R]
    elseif haskey(hint.t, (jname, iname, -R))
      t[tidx] = conj(hint.t[jname, iname, -R])
    end
  end

  Γ, Δ = HFB.computetargetfields(computer, ρ, t)
  return HFBSolution(ρ, t, Γ, Δ)
end


function newhfbhint(computer::HFBComputer{O}, sol::HFBSolution) where {O}
  ρ = Dict{Tuple{O, O, Vector{Int64}}, Complex128}()
  t = Dict{Tuple{O, O, Vector{Int64}}, Complex128}()
  uc = computer.unitcell

  for (ρidx, (i, j, rij)) in enumerate(computer.ρ_registry)
    iname, icoord = getorbital(uc, i)
    jname, jcoord = getorbital(uc, j)
    Ri = whichunitcell(uc, iname, fract2carte(uc, icoord))
    Rj = whichunitcell(uc, jname, fract2carte(uc, icoord) + rij)
    @assert(!haskey(ρ, (iname, jname, Rj-Ri)))
    ρ[iname,jname,Rj-Ri] = sol.ρ[ρidx]
  end

  for (tidx, (i, j, rij)) in enumerate(computer.t_registry)
    iname, icoord = getorbital(uc, i)
    jname, jcoord = getorbital(uc, j)
    Ri = whichunitcell(uc, iname, fract2carte(uc, icoord))
    Rj = whichunitcell(uc, jname, fract2carte(uc, icoord) + rij)
    @assert(!haskey(t, (iname, jname, Rj-Ri)))
    t[iname,jname,Rj-Ri] = sol.t[tidx]
  end

  return HFBHint{O}(ρ, t)
end


function fixhfbsolution(computer::HFBComputer{O},
                        sol ::HFBSolution) where {O}
  sol.Γ[:], sol.Δ[:] = HFB.computetargetfields(computer, sol.ρ, sol.t)
end


function randomize!(computer::HFBComputer{O},
                    sol ::HFBSolution,
                    amplitude::C=1.0) where {O,C<:Number}
  sol.ρ[:] = (rand(Float64, length(sol.ρ)) .* 2 .- 1) * amplitude
  sol.t[:] = (rand(Complex128, length(sol.t)) .* 2 .- 1) * amplitude
  sol.Γ[:], sol.Δ[:] = HFB.computetargetfields(computer, sol.ρ, sol.t)
end
