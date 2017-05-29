module HartreeFockBogoliubovSolver

using ..Lattice
import ..Embed
import ..Spec


type MeanField
  target ::Tuple{Int64, Int64}
  source ::Tuple{Int64, Int64}
  amplitude ::Real
  targetdisplacement ::Vector{Float64}
  sourcedisplacement ::Vector{Float64}
end


# Nambu Space
type HartreeFockBogoliubovModel
  hoppings #:: function which returns a matrix.Sparse or Dense?
  particle_hole_interactions
  particle_particle_interactions
end


function addinteraction!(hfbmodel ::HartreeFockBogoliubovModel,
                         embint ::Embed.InteractionDiagonal)
  V = embint.amplitude
  (i, j)  = (embint.i,  embint.j )
  (ri,rj) = (embint.ri, embint.rj)

  Γ = [
    MeanField((i,i), (j,j),      V, ri - ri, rj - rj),
    MeanField((j,i), (i,j),     -V, rj - ri, ri - rj),
    MeanField((i,j), (j,i),     -V, ri - rj, rj - ri),
    MeanField((j,j), (i,i),      V, rj - rj, ri - ri),
  ]
  Δ = [
    MeanField((i,j), (i,j),  0.5*V, ri - rj, ri - rj),
    MeanField((i,j), (j,i), -0.5*V, ri - rj, rj - ri),
    MeanField((j,i), (i,j), -0.5*V, rj - ri, ri - rj),
    MeanField((j,i), (j,i),  0.5*V, rj - ri, rj - ri),
  ]
  append!(hfbmodel.particle_hole_interactions, Γ)
  append!(hfbmodel.particle_particle_interactions, Δ)
end



function addinteraction!(hfbmodel ::HartreeFockBogoliubovModel,
                         embint ::Embed.InteractionOffdiagonal)
  error("Unimplemented")
  V = embint.amplitude
  Vc = conjugate(V)
  (i, j, k, l)  = (embint.i,  embint.j,  embint.k,  embint.l )
  (ri,rj,rk,rl) = (embint.ri, embint.rj, embint.rk, embint.rl)

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
  append!(hfbmodel.particle_hole_interactions, Γs)
  append!(hfbmodel.particle_particle_interactions, Δs)
end



function makecollector(hfbmodel ::HartreeFockBogoliubovModel)
  # TODO(kmlee, 2017/05/29) Policy: unit cell, or cartesian?
  ρ_collect = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in hfbmodel.particle_hole_interaction
    (k, l) = Gamma.source
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCkl = fract2carte(hfbmodel.unitcell.latticevectors, rFkl)

    if haskey(ρ_collect, (k, l, rFkl.w))
      prev_func = ρ_collect[(k, l, rFkl.w)]
      ρ_collect[(k, l, rFkl.w)] = (momentum::Vector{Float64}, ρ) -> begin
        return (ρ(k,l) * exp(-1im * dot(rCkl, momentum)) + prev_func(momentum, ρ))
      end
    else
      ρ_collect[(k, l, rFkl.w)] = (momentum::Vector{Float64}, ρ) -> begin
        return (ρ(k,l) * exp(-1im * dot(rCkl, momentum)))
      end
    end
  end

  t_collect = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in hfbmodel.particle_particle_interaction
    (k, l) = Gamma.source
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCkl = fract2carte(hfbmodel.unitcell.latticevectors, rFkl)

    if haskey(t_collect, (k, l, rFkl.w))
      prev_func = t_collect[(k, l, rFkl.w)]
      t_collect[(k, l, rFkl.w)] = (momentum::Vector{Float64}, t) -> begin
        return (t(k,l) * exp(-1im * dot(rCkl, momentum)) + prev_func(momentum, t))
      end
    else
      t_collect[(k, l, rFkl.w)] = (momentum::Vector{Float64}, t) -> begin
        return (t(k,l) * exp(-1im * dot(rCkl, momentum)))
      end
    end
  end
  return (ρ_collect, t_collect)
end



function collect(unitcell ::UnitCell, ρ_collect, t_collect, momentum, ρ, t)
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



function makedeployer(hfbmodel ::HartreeFockBogoliubovModel)
  ρ_deploy = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in hfbmodel.particle_hole_interaction
    (i, j) = Gamma.target
    (k, l) = Gamma.source
    V = Gamma.amplitude
    rFij::FractCoord = Gamma.targetdisplacement
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCij = fract2carte(hfbmodel.unitcell.latticevectors, rFij)

    if haskey(ρ_deploy, (i, j, rFij.w))
      prev_func = ρ_deploy[(i, j, rFij.w)]
      ρ_deploy[(i, j, rFij.w)] = (momentum, ρ_value) -> begin
        return (V * ρ_value[(k, l, rFkl.w)] * exp(1im * dot(rCij, momentum))
                + prev_func(momentum, ρ_value))
      end
    else
      ρ_deploy[(i, j, rFij.w)] = (momentum, ρ_value) -> begin
        return (V * ρ_value[(k, l, rFkl.w)] * exp(1im * dot(rCij, momentum)))
      end
    end
  end

  t_deploy = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in hfbmodel.particle_particle_interaction
    (i, j) = Gamma.target
    (k, l) = Gamma.source
    V = Gamma.amplitude
    rFij::FractCoord = Gamma.targetdisplacement
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCij = fract2carte(hfbmodel.unitcell.latticevectors, rFij)

    if haskey(t_deploy, (i, j, rFij.w))
      prev_func = t_deploy[(i, j, rFij.w)]
      t_deploy[(i, j, rFij.w)] = (momentum, t_value) -> begin
        return (V * t_value[(k, l, rFkl.w)] * exp(1im * dot(rCij, momentum))
                + prev_func(momentum, t_value))
      end
    else
      t_deploy[(i, j, rFij.w)] = (momentum, t_value) -> begin
        return (V * t_value[(k, l, rFkl.w)] * exp(1im * dot(rCij, momentum)))
      end
    end
  end
  return (ρ_deploy, t_deploy)
end


function deploy(unitcell ::UnitCell, ρ_deploy, t_deploy, momentum, ρ_value, ρ_value)
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


# Nambu Space
type HartreeFockBogoliubovModelRealization
  hoppings #:: function which returns a matrix.Sparse or Dense?
  particle_hole_interactions
  particle_particle_interactions
end









end