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
       fixhfbsolution,
       iscompatible

export makehoppingmatrix,
       makeGammamatrix,
       makeDeltamatrix

export freeze

using HartreeFockBogoliubov

"""
`CollectRow` is holds info on how to compute ρ or t.
Its elements are:
1. Is diagonal? (only for rho)
2. row orbital
3. col orbital
4. displacement r(col) - r(row)
"""
const CollectRow = Tuple{Bool, Int64, Int64, Vector{Float64}}


"""
`DeployRow` is holds info on how to compute Γ or Δ.
Its elements are:
1. Is diagonal
2. row orbital
3. col orbital
4. displacement r(col) - r(row)
5. list of sources, each of which is a tuple of
  1. index of ρ or t from which to compute this Γ or Δ.
  2. amplitude (coefficient to multiply to ρ or t)
  3. boolean indicating whether
     (1) conjugation is needed (for ρ/Γ) or
     (2) minus sign is needed (for t/Δ).
"""
const DeployRow = Tuple{Bool, Int64, Int64, Vector{Float64}, Vector{Tuple{Int64, Complex128, Bool}}}


"""
`HFBConmputer` is a type holding the ρ, t and Γ, Δ of a Hartree-Fock-Bogoliubov
Hamiltonian.
"""
mutable struct HFBComputer{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{Spec.Hopping}
  temperature ::Float64

  fermi ::Function
  ρ_registry ::Vector{CollectRow}
  t_registry ::Vector{CollectRow}
  Γ_registry ::Vector{DeployRow}
  Δ_registry ::Vector{DeployRow}
end


"""
"""
function HFBComputer(unitcell::UnitCell{O},
                     hoppings::AbstractVector{Spec.Hopping},
                     temperature::Real) where {O}
 fermi = fermidirac(temperature)
 return HFBComputer{O}(unitcell, hoppings, temperature, fermi, [], [], [], [])
end


"""
"""
function makeparticleholeregistry(ham::HFBHamiltonian{T}) where {T}
  function getdistance(i ::Integer, j ::Integer, Rij::AbstractVector{<:Integer})
    ri, rj = getorbitalcoord(ham.unitcell, i), getorbitalcoord(ham.unitcell, j)
    rj = rj + Rij
    ri, rj = fract2carte(ham.unitcell, ri), fract2carte(ham.unitcell, rj)
    return rj - ri
  end

  collect_reg = Dict()
  deploy_reg = Dict()

  let
    R = zeros(Int64, dimension(ham.unitcell))
    r = zeros(Float64, dimension(ham.unitcell))
    for (i, orb) in enumerate(ham.unitcell.orbitals)
      collect_reg[i, i, zeros(R)] = (length(collect_reg)+1, (true, i, i, zeros(r)))
    end
    for (i, orb) in enumerate(ham.unitcell.orbitals)
      deploy_reg[i, i, zeros(R)] = (length(deploy_reg)+1, (true, i, i, zeros(r), []))
    end
  end

  for hopmf in ham.particle_hole_interactions
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    rij = getdistance(i, j, Rij)
    rkl = getdistance(k, l, Rkl)
    if !haskey(deploy_reg, (i,j, Rij))
      isdiag = (i == j) && iszero(Rij)
      deploy_reg[i,j,Rij] = (length(deploy_reg)+1, (isdiag, i, j, rij, []))
    end
    if !haskey(collect_reg, (k,l, Rkl))
      isdiag = (k == l) && iszero(Rkl)
      collect_reg[k,l,Rkl] = (length(collect_reg)+1, (isdiag, k, l, rkl))
    end
  end

  # Add sources to targets
  for hopmf in ham.particle_hole_interactions
    v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    srcidx = collect_reg[k,l,Rkl][1]
    star = hopmf.star
    if abs(v) > eps(Float64)
      push!( deploy_reg[i,j,Rij][2][5], (srcidx, v, star) )
    end
  end

  let indices = [idx for (key, (idx, val)) in collect_reg]
    @assert(length(unique(indices)) == length(indices))
  end
  let indices = [idx for (key, (idx, val)) in deploy_reg]
    @assert(length(unique(indices)) == length(indices))
  end

  ρ_registry = sort([(idx, val) for (key, (idx, val)) in collect_reg], by=(x) -> x[1])
  Γ_registry = sort([(idx, val) for (key, (idx, val)) in deploy_reg ], by=(x) -> x[1])
  ρ_registry = [val for (idx, val) in ρ_registry]
  Γ_registry = [val for (idx, val) in Γ_registry]

  return (ρ_registry, Γ_registry)
end


"""
"""
function makeparticleparticleregistry(ham::HFBHamiltonian{T}) where {T}
  function getdistance(i ::Integer, j ::Integer, Rij::AbstractVector{<:Integer})
    ri, rj = getorbitalcoord(ham.unitcell, i), getorbitalcoord(ham.unitcell, j)
    rj = rj + Rij
    ri, rj = fract2carte(ham.unitcell, ri), fract2carte(ham.unitcell, rj)
    return rj - ri
  end

  collect_reg = Dict()
  deploy_reg = Dict()

  for hopmf in ham.particle_particle_interactions
    #v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    rij = getdistance(i, j, Rij)
    rkl = getdistance(k, l, Rkl)
    @assert( !(i==j && iszero(Rij)) )
    @assert( !(k==l && iszero(Rkl)) )
    if !haskey(collect_reg, (k,l,Rkl))
      collect_reg[k,l,Rkl] = (length(collect_reg)+1, (false, k, l, rkl))
    end
    if !haskey(deploy_reg, (i,j,Rij))
      deploy_reg[i,j,Rij] = (length(deploy_reg)+1, (false, i, j, rij, []))
    end
  end

  for hopmf in ham.particle_particle_interactions
    v = hopmf.amplitude
    (i,j,Rij) = hopmf.target
    (k,l,Rkl) = hopmf.source
    srcidx = collect_reg[k,l,Rkl][1]
    neg = hopmf.negate
    if abs(v) > eps(Float64)
      push!( deploy_reg[i,j,Rij][2][5], (srcidx, v, neg) )
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
  return (t_registry, Δ_registry)
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

  (ρ_registry, Γ_registry) = makeparticleholeregistry(ham)
  (t_registry, Δ_registry) = makeparticleparticleregistry(ham)

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

function iscompatible(s1 ::HFBSolution, s2::HFBSolution)
  return (length(s1.ρ) == length(s2.ρ) &&
          length(s1.t) == length(s2.t) &&
          length(s1.Γ) == length(s2.Γ) &&
          length(s1.Δ) == length(s2.Δ))
end


import Base: copy

copy(x::HFBSolution) = HFBSolution(copy(x.ρ), copy(x.t), copy(x.Γ), copy(x.Δ))
copy{T}(x::HFBHint{T}) = HFBHint{T}(copy(x.ρ), copy(x.t))

import Base: +, -, *, /, \

function +(sol1::HFBSolution)
  HFBSolution(sol1.ρ, sol1.t, sol1.Γ, sol1.Δ)
end

function -(sol1::HFBSolution)
  HFBSolution(-sol1.ρ, -sol1.t, -sol1.Γ, -sol1.Δ)
end

function +(sol1::HFBSolution, sol2::HFBSolution)
  HFBSolution(sol1.ρ + sol2.ρ, sol1.t + sol2.t,
              sol1.Γ + sol2.Γ, sol1.Δ + sol2.Δ)
end

function -(sol1::HFBSolution, sol2::HFBSolution)
  HFBSolution(sol1.ρ - sol2.ρ, sol1.t - sol2.t,
              sol1.Γ - sol2.Γ, sol1.Δ - sol2.Δ)
end

function *(sol1::HFBSolution, value :: T) where {T <: Number}
  HFBSolution(sol1.ρ * value, sol1.t * value,
              sol1.Γ * value, sol1.Δ * value)
end

function *(value :: T, sol1::HFBSolution) where {T <: Number}
  HFBSolution(value * sol1.ρ, value * sol1.t,
              value * sol1.Γ, value * sol1.Δ)
end

function /(sol1::HFBSolution, value :: T) where {T <: Number}
  HFBSolution(sol1.ρ / value, sol1.t / value,
              sol1.Γ / value, sol1.Δ / value)
end

function \(value :: T, sol1::HFBSolution) where {T <: Number}
  HFBSolution(value \ sol1.ρ, value \ sol1.t,
              value \ sol1.Γ, value \ sol1.Δ)
end


"""
func : (idx, i, j, r) -> val
"""
function makesourcefields(funcρ ::Function,
                          funct ::Function,
                          computer ::HFBComputer{O}) where {O}
  ρs = zeros(Complex128, length(computer.ρ_registry))
  ts = zeros(Complex128, length(computer.t_registry))

  for (idx, (isdiag, i, j, r)) in enumerate(computer.ρ_registry)
    v = funcρ(idx, i, j, r)
    #if i==j && all((x)->x==0, r)
    if isdiag
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


"""
func : (idx, i, j, r) -> 0
"""
function makesourcefields(computer ::HFBComputer{O}) where {O}
  ρs = zeros(Complex128, length(computer.ρ_registry))
  ts = zeros(Complex128, length(computer.t_registry))
  return (ρs, ts)
end

"""
Compute Γ and Δ from ρ and t.
"""
function computetargetfields(computer ::HFBComputer{O},
                             ρs ::AbstractVector{<:Number},
                             ts ::AbstractVector{<:Number}) where {O}
  Γs = zeros(Complex128, length(computer.Γ_registry))
  Δs = zeros(Complex128, length(computer.Δ_registry))
  for (tgtidx, (_, i, j, r, srcs)) in enumerate(computer.Γ_registry)
    value = 0.0 + 0.00im
    for (srcidx, amplitude, star) in srcs
      value += amplitude * (star ? conj(ρs[srcidx]) : ρs[srcidx])
    end
    Γs[tgtidx] = value
  end

  for (tgtidx, (_, i, j, r, srcs)) in enumerate(computer.Δ_registry)
    value = 0.0 + 0.00im
    for (srcidx, amplitude, neg) in srcs
      value += amplitude * (neg ? -ts[srcidx] : ts[srcidx])
    end
    Δs[tgtidx] = value
  end
  return (Γs, Δs)
end

"""
Return a generator of hopping matrix (which is a function of momentum)
"""
function makehoppingmatrix(computer::HFBComputer)
  norb = numorbital(computer.unitcell)
  hk = Generator.generatefast(computer.unitcell, computer.hoppings)
  function(k ::AbstractVector{Float64})
    out = zeros(Complex128, (norb, norb))
    hk( k, out )
    return out
  end
end

"""
Return a generator of Γ matrix (which is a function of momentum)
"""
function makeGammamatrix(computer::HFBComputer,
                         Γs ::AbstractVector{<:Number})
  norb = numorbital(computer.unitcell)
  function(k ::AbstractVector{Float64})
    out = zeros(Complex128, (norb, norb))
    for (idx, Γ) in enumerate(Γs)
      (isdiag, i, j, r, _) = computer.Γ_registry[idx]
      if isdiag
        @assert(i==j && all(x -> isapprox(x, 0.0), r))
        @assert(isapprox(imag(Γ), 0.0))
        out[i,j] += real(Γ)
      else
        @assert(!(i==j && all(x -> isapprox(x, 0.0), r)))
        phase = cis(dot( k, r))
        out[i,j] += Γ * phase
        out[j,i] += conj(Γ * phase)
      end
    end
    return out
  end
end


"""
Return a generator of Δ matrix (which is a function of momentum)
"""
function makeDeltamatrix(computer::HFBComputer,
                         Δs ::AbstractVector{<:Number})
  norb = numorbital(computer.unitcell)
  function(k ::AbstractVector{Float64})
    out = zeros(Complex128, (norb, norb))
    # 3/3. Delta
    for (idx, Δ) in enumerate(Δs)
      (isdiag, i, j, r, _) = computer.Δ_registry[idx]
      @assert( !isdiag )
      @assert( !(i==j && all(x -> isapprox(x, 0.0), r)) )
      phase = cis(dot(k, r))
      out[i,j] += Δ * phase
      out[j,i] -= Δ * conj(phase)
    end
    return out
  end
end


"""
"""
function makehamiltonian(computer ::HFBComputer,
                         Γs ::AbstractVector{<:Number},
                         Δs ::AbstractVector{<:Number})
  norb = numorbital(computer.unitcell)
  hk = Generator.generatefast(computer.unitcell, computer.hoppings)

  function(k ::AbstractVector{Float64})
    out = zeros(Complex128, (norb, 2, norb, 2))
    # 1/3. non-interacting kinetic part
    hk( k, view(out, :,1,:,1))
    hk(-k, view(out, :,2,:,2))
    # 2/3. Gamma
    for (idx, Γ) in enumerate(Γs)
      (isdiag, i, j, r, _) = computer.Γ_registry[idx]
      if isdiag
        @assert(i==j && all(x -> isapprox(x, 0.0), r))
        @assert(isapprox(imag(Γ), 0.0))
        out[i,1,j,1] += real(Γ)
        out[i,2,j,2] += real(Γ)
      else
        @assert(!(i==j && all(x -> isapprox(x, 0.0), r)))
        phase = cis(dot( k, r))
        out[i,1,j,1] += Γ * phase
        out[i,2,j,2] += Γ * conj(phase)
        out[j,1,i,1] += conj(Γ * phase)
        out[j,2,i,2] += conj(Γ) * phase
      end
    end
    out[:,2,:,2] = -transpose(out[:,2,:,2])

    # 3/3. Delta
    for (idx, Δ) in enumerate(Δs)
      (isdiag, i, j, r, _) = computer.Δ_registry[idx]
      @assert( !isdiag )
      @assert( !(i==j && all(x -> isapprox(x, 0.0), r)) )
      phase = cis(dot(k, r))
      val1 = Δ * phase
      val2 = Δ * conj(phase)
      out[i,1,j,2] += val1
      out[j,1,i,2] -= val2
      out[j,2,i,1] += conj(val1)
      out[i,2,j,1] -= conj(val2)
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

  function(k::AbstractVector{<:Real},
           eigenvalues ::AbstractVector{<:Real},
           eigenvectors ::AbstractMatrix{<:Number},
           ρout ::AbstractVector{Complex128},
           tout ::AbstractVector{Complex128})
    @assert(length(eigenvalues) == 2*norb)
    @assert(size(eigenvectors) == (2*norb, 2*norb))

    f = [fermi(e) for e in eigenvalues]
    ψ = reshape(eigenvectors, (norb, 2, norb*2))
    u = ψ[:, 1, :]
    v = ψ[:, 2, :]

    ρfunc(i ::Int, j ::Int) = sum(f .* u[i, :] .* conj(u[j, :]))
    tfunc(i ::Int, j ::Int) = sum(f .* u[i, :] .* conj(v[j, :]))

    for (idx, (isdiag, i, j, r)) in enumerate(ρ_registry)
      if isdiag
        @assert((i==j && all(x -> isapprox(x, 0.0), r)))
        ρout[idx] += real( ρfunc(i, i) )
      else
        @assert(!(i==j && all(x -> isapprox(x, 0.0), r)))
        ρout[idx] += ρfunc(i, j) * cis(-dot(k, r))
      end
    end
    for (idx, (isdiag, i, j, r)) in enumerate(t_registry)
      @assert(!isdiag)
      @assert(!(i==j && all(x -> isapprox(x, 0.0), r)))
      tout[idx] += tfunc(i, j) * cis(-dot(k, r))
    end
  end
end


"""
  Return a zero solution
"""
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


"""
"""
function newhfbhint(computer::HFBComputer{O}, sol::HFBSolution) where {O}
  ρ = Dict{Tuple{O, O, Vector{Int64}}, Complex128}()
  t = Dict{Tuple{O, O, Vector{Int64}}, Complex128}()
  uc = computer.unitcell

  for (ρidx, (_, i, j, rij)) in enumerate(computer.ρ_registry)
    iname, icoord = getorbital(uc, i)
    jname, jcoord = getorbital(uc, j)
    Ri = whichunitcell(uc, iname, fract2carte(uc, icoord))
    Rj = whichunitcell(uc, jname, fract2carte(uc, icoord) + rij)
    @assert(!haskey(ρ, (iname, jname, Rj-Ri)))
    ρ[iname,jname,Rj-Ri] = sol.ρ[ρidx]
  end

  for (tidx, (_, i, j, rij)) in enumerate(computer.t_registry)
    iname, icoord = getorbital(uc, i)
    jname, jcoord = getorbital(uc, j)
    Ri = whichunitcell(uc, iname, fract2carte(uc, icoord))
    Rj = whichunitcell(uc, jname, fract2carte(uc, icoord) + rij)
    @assert(!haskey(t, (iname, jname, Rj-Ri)))
    t[iname,jname,Rj-Ri] = sol.t[tidx]
  end

  return HFBHint{O}(ρ, t)
end


"""
    Recompute Γ and Δ from ρ and t in a `HFBSolution`
"""
function fixhfbsolution(computer::HFBComputer{O},
                        sol ::HFBSolution) where {O}
  sol.Γ[:], sol.Δ[:] = HFB.computetargetfields(computer, sol.ρ, sol.t)
end


"""
    Randomize a solution
"""
function randomize!(computer::HFBComputer{O},
                    sol ::HFBSolution) where {O}
  sol.ρ[:] = (rand(Float64, length(sol.ρ)))
  sol.t[:] = (rand(Complex128, length(sol.t)))
  sol.Γ[:], sol.Δ[:] = HFB.computetargetfields(computer, sol.ρ, sol.t)
end


"""
"""
function nambufy(uc::UnitCell{O}, hop::HoppingDiagonal{R}) where {O, R<:Real}
  norb = numorbital(uc)
  hop1 = HoppingDiagonal{R}( hop.amplitude, hop.i, hop.Ri)
  hop2 = HoppingDiagonal{R}(-hop.amplitude, hop.i+norb, hop.Ri)
  return (hop1, hop2)
end


"""
"""
function nambufy(uc::UnitCell{O}, hop::HoppingOffdiagonal{C}) where {O, C<:Number}
  norb = numorbital(uc)
  hop1 = HoppingOffdiagonal{C}(      hop.amplitude,  hop.i     , hop.j     , hop.Ri, hop.Rj)
  hop2 = HoppingOffdiagonal{C}(-conj(hop.amplitude), hop.i+norb, hop.j+norb, hop.Ri, hop.Rj)
  return (hop1, hop2)
end



"""
"""
function nambufy(uc::UnitCell{O}, hop::HoppingDiagonal{R}) where {O, R<:Real}
  norb = numorbital(uc)
  hop1 = HoppingDiagonal{R}( hop.amplitude, hop.i, hop.Ri)
  hop2 = HoppingDiagonal{R}(-hop.amplitude, hop.i+norb, hop.Ri)
  return (hop1, hop2)
end


"""
"""
function nambufy(uc::UnitCell{O}, hop::HoppingOffdiagonal{C}) where {O, C<:Number}
  norb = numorbital(uc)
  hop1 = HoppingOffdiagonal{C}(      hop.amplitude,  hop.i     , hop.j     , hop.Ri, hop.Rj)
  hop2 = HoppingOffdiagonal{C}(-conj(hop.amplitude), hop.i+norb, hop.j+norb, hop.Ri, hop.Rj)
  return (hop1, hop2)
end


"""
    freeze

    Create a hopping hamiltonian out of HFB Hamiltonian and Solution

    Order: 
    electrons......, HOLE......
"""
function freeze(computer ::HFBComputer{O},
                Γs ::AbstractVector{<:Number},
                Δs ::AbstractVector{<:Number}) where {O}
  norb = numorbital(computer.unitcell)
  unitcell = computer.unitcell

  NewOrbitalType = (O <: Tuple) ? Tuple{O.parameters..., Symbol} : Tuple{O, Symbol}
  nambuunitcell = Lattice.newunitcell(unitcell.latticevectors; OrbitalType=NewOrbitalType)
  if O <: Tuple
    for (orb, fc) in unitcell.orbitals
      addorbital!(nambuunitcell, (orb..., :PARTICLE), fc)
    end
    for (orb, fc) in unitcell.orbitals
      addorbital!(nambuunitcell, (orb..., :HOLE), fc)
    end
  else
    for (orb, fc) in unitcell.orbitals
      addorbital!(nambuunitcell, (orb, :PARTICLE), fc)
    end
    for (orb, fc) in unitcell.orbitals
      addorbital!(nambuunitcell, (orb, :HOLE), fc)
    end
  end
  hoppings = Hopping[]
  for hop in computer.hoppings
    (hop1, hop2) = nambufy(unitcell, hop)
    push!(hoppings, hop1)
    push!(hoppings, hop2)
  end

  # 2/3. Gamma
  for (idx, Γ) in enumerate(Γs)
    (isdiag, i, j, r, _) = computer.Γ_registry[idx]
    if isdiag
      @assert(i==j && all(x -> isapprox(x, 0.0), r))
      @assert(isapprox(imag(Γ), 0.0))

      iname, icoord = getorbital(unitcell, i)
      Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))

      push!(hoppings, HoppingDiagonal{Float64}( real(Γ), i, Ri))
      push!(hoppings, HoppingDiagonal{Float64}(-real(Γ), i+norb, Ri))
    else
      @assert(!(i==j && all(x -> isapprox(x, 0.0), r)))
      
      iname, icoord = getorbital(unitcell, i)
      jname, jcoord = getorbital(unitcell, j)
      Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
      Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)

      push!(hoppings, HoppingOffdiagonal{Complex128}( Γ      , i, j, Ri, Rj))
      push!(hoppings, HoppingOffdiagonal{Complex128}(-conj(Γ), i+norb, j+norb, Ri, Rj))
    end
  end

  # 3/3. Delta
  for (idx, Δ) in enumerate(Δs)
    (isdiag, i, j, r, _) = computer.Δ_registry[idx]

    @assert( !isdiag )
    @assert( !(i==j && all(x -> isapprox(x, 0.0), r)) )
    
    iname, icoord = getorbital(unitcell, i)
    jname, jcoord = getorbital(unitcell, j)
    Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
    Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)

    push!(hoppings, HoppingOffdiagonal{Complex128}( Δ , i, j + norb, Ri, Rj))
    push!(hoppings, HoppingOffdiagonal{Complex128}(-Δ , j, i + norb, Rj, Ri))
  end
  return (nambuunitcell, hoppings)
end
