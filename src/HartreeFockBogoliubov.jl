module HartreeFockBogoliubov

export HartreeFockBogoliubovModel
export addinteraction
export MeanField

type MeanField
  target ::Tuple{Int64, Int64}
  source ::Tuple{Int64, Int64}
  amplitude ::Real
  displacement ::Vector{Float64}
end

# Nambu Space
type HartreeFockBogoliubovModel
  hoppings #:: function which returns a matrix.Sparse or Dense?
  particle_hole_interactions
  particle_particle_interactions
end

"""
```math
  H =  \sum_{ij} T_{ij} c_{i}^{\dagger} c_{j}
     + \frac{1}{4} \sum_{ijkl} V_{ijkl} c_{i}^{\dagger} c_{j}^{\dagger} c_{l} c_{k}
````

```math
  \Gamma_{ij} = \sum_{kl} V_{ikjl} \rho_{kl}
```

```math
  \Delta_{ij} = 1/2 \sum_{kl} V_{ijkl} t_{kl}
``

  <---->
^ i    j ^
|        |
|        |
v l    k v
  <---->
"""
function addinteraction(hfbmodel ::HartreeFockBogoliubovModel,
                        amplitude::Real,
                        i::Integer, j::Integer,
                        k::Integer, l::Integer,
                        )
  Γ = [
    MeanField((i,k), (j,l), amplitude, displacement)
  ]
  Γ = [
    ((i,k), (j,l),  amplitude),
    ((j,k), (i,l), -amplitude),
    ((i,l), (j,k), -amplitude),
    ((j,l), (i,k),  amplitude),
  ]
  Δ = [
    ((i,j), (k,l),  0.5*amplitude),
    ((i,j), (l,k), -0.5*amplitude),    
    ((j,i), (k,l), -0.5*amplitude),
    ((j,i), (l,k),  0.5*amplitude),    
  ]
  append!(hfbmodel.particle_hole_interactions, Γ)
  append!(hfbmodel.particle_particle_interactions, Δ)
end


function make_hopping_generator_maker(hoppings :: Vector)
end

# Nambu Space
type HartreeFockBogoliubovModelRealization
  hoppings #:: function which returns a matrix.Sparse or Dense?
  particle_hole_interactions
  particle_particle_interactions
end


end