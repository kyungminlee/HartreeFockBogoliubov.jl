type HartreeFockBogoliubovSolver
  unitcell ::UnitCell
  hoppings ::Vector{Embed.Hopping}
  ρ_collectors ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
  t_collectors ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
  ρ_deployers  ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
  t_deployers  ::Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}
end

function HartreeFockBogoliubovSolver(model ::Hamiltonian,
                                     unitcell::UnitCell)
  lv ::Array{Float64, 2} = unitcell.latticevectors
  ph = model.particle_hole_interactions
  pp = model.particle_particle_interactions

  ρ_collectors = makecollectors(ph)
  t_collectors = makecollectors(pp)
  ρ_deployers  = makedeployers(ph, lv)
  t_deployers  = makedeployers(pp, lv)
  return HartreeFockBogoliubovSolver(model.hoppings,
                                     ρ_collectors, t_collectors,
                                     ρ_deployers, t_deployers)
end


"""
"""
function makecollectors(interactions ::Vector{Embed.Interaction})
  collectors = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in interactions
    (k, l) = Gamma.source
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCkl::CarteCoord = fract2carte(hfbmodel.unitcell.latticevectors, rFkl)

    if haskey(collecters, (k, l, rFkl.w))
      prev_func = collectors[(k, l, rFkl.w)]
      collectors[(k, l, rFkl.w)] = (momentum::Vector{Float64}, green) -> begin
        return (green(k,l) * exp(-1im * dot(rCkl, momentum))
                + prev_func(momentum, green))
      end
    else
      collectors[(k, l, rFkl.w)] = (momentum::Vector{Float64}, green) -> begin
        return (green(k,l) * exp(-1im * dot(rCkl, momentum)))
      end
    end
  end
  return collectors
end


"""
"""
function makedeployers(interactions ::Vector{Embed.Interaction},
                       latticevectors ::Array{Float64, 2})
  deployers = Dict{Tuple{Int64, Int64, Vector{Int64}}, Any}()
  for Gamma in interactions
    (i, j) = Gamma.target
    (k, l) = Gamma.source
    V = Gamma.amplitude
    rFij::FractCoord = Gamma.targetdisplacement
    rFkl::FractCoord = Gamma.sourcedisplacement
    rCij::CarteCoord = fract2carte(latticevectors, rFij)

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
