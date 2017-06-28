export HFBHamiltonian
export addinteraction!


"""
"""
struct HoppingMeanField{C<:Number}
  amplitude ::C
  target ::Tuple{Int64, Int64, Vector{Int64}}
  source ::Tuple{Int64, Int64, Vector{Int64}}
  star ::Bool

  function HoppingMeanField{C}(
          amplitude ::C,
          target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
          source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
          star ::Bool
        ) where {C<:Number}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    @assert(i <= j)
    @assert(k <= l)
    new{C}(amplitude, target, source, star)
  end

  function HoppingMeanField(
          amplitude ::C,
          target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
          source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}}
        ) where {C<:Number}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    star = false
    if i > j
      (i,j) = (j,i)
      Rij = -Rij
      amplitude = conj(amplitude)
      star = !star
    end
    if k > l
      (k,l) = (l,k)
      Rkl = -Rkl
      star = !star
    end
    new{C}(amplitude, (i,j,Rij), (k,l,Rkl), star)
  end
end


"""
"""
struct PairingMeanField{C<:Number}
  amplitude ::C
  target ::Tuple{Int64, Int64, Vector{Int64}}
  source ::Tuple{Int64, Int64, Vector{Int64}}
  negate ::Bool

  function PairingMeanField{C}(
          amplitude ::C,
          target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
          source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
          negate ::Bool
        ) where {C<:Number}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    @assert(i <= j)
    @assert(k <= l)
    new{C}(amplitude, target, source, negate)
  end

  function PairingMeanField(
          amplitude ::C,
          target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
          source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}}
        ) where {C<:Number}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    negate = false
    if i > j
      (i,j) = (j,i)
      Rij = -Rij
      amplitude = -amplitude
      #TODO HERE!!! flip amplitude or not flip amplitude
      #negate = !negate
    end
    if k > l
      (k,l) = (l,k)
      Rkl = -Rkl
      amplitude = -amplitude
      #negate = !negate
    end
    new{C}(amplitude, (i,j,Rij), (k,l,Rkl), negate)
  end
end

"""
"""
mutable struct HFBHamiltonian{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{Hopping}
  particle_hole_interactions ::Vector{HoppingMeanField}
  particle_particle_interactions ::Vector{PairingMeanField}
end


"""
"""
function HFBHamiltonian(
        unitcell ::UnitCell{O},
        hoppings ::AbstractVector{Hopping},
        particle_hole_interactions ::AbstractVector{HoppingMeanField},
        particle_particle_interactions ::AbstractVector{PairingMeanField}
      ) where {O}
  return HFBHamiltonian{O}(unitcell,
                           hoppings,
                           particle_hole_interactions,
                           particle_particle_interactions)
end


"""
"""
function HFBHamiltonian(full_hamiltonian ::FullHamiltonian{O}) where {O}
  unitcell = full_hamiltonian.unitcell
  hoppings = full_hamiltonian.hoppings
  model = HFBHamiltonian{O}(unitcell, hoppings, [], [])
  for interaction in full_hamiltonian.interactions
    addinteraction!(model, interaction)
  end
  return model
end


"""
Add diagonal interaction
"""
function addinteraction!(model ::HFBHamiltonian{O},
                         specint ::InteractionDiagonal{R}) where {O,R<:Real}
  v = specint.amplitude
  (i, j) = (specint.i,  specint.j )
  (Ri, Rj) = (specint.Ri, specint.Rj)

  Γs = [
    HoppingMeanField( v, (i, i, Ri-Ri), (j, j, Rj-Rj)),
    HoppingMeanField( v, (j, j, Rj-Rj), (i, i, Ri-Ri)),
    HoppingMeanField(-v, (i, j, Rj-Ri), (i, j, Rj-Ri)),
  ]
  Δs = [
    PairingMeanField( v, (i, j, Rj-Ri), (i, j, Rj-Ri)),
  ]
  append!(model.particle_hole_interactions, Γs)
  append!(model.particle_particle_interactions, Δs)
end


"""
Add diagonal interaction
"""
function addinteraction!(model ::HFBHamiltonian{O},
                         specint ::InteractionOffdiagonal{C}) where {O,C<:Number}
  v = specint.amplitude
  (i, j) = (specint.i,  specint.j)
  (k, l) = (specint.k,  specint.l)
  (Ri, Rj) = (specint.Ri, specint.Rj)
  (Rk, Rl) = (specint.Rk, specint.Rl)

  Γs = [
    HoppingMeanField( v, (i, k, Rk-Ri), (l, j, Rj-Rl)),
    HoppingMeanField(-v, (i, l, Rl-Ri), (k, j, Rj-Rk)),
    HoppingMeanField(-v, (j, k, Rk-Rj), (l, i, Ri-Rl)),
    HoppingMeanField( v, (j, l, Rl-Rj), (k, i, Ri-Rk)),
  ]
  Δs = [
    PairingMeanField( v, (i, j, Rj-Ri), (k, l, Rl-Rk)),
    PairingMeanField( conj(v), (k, l, Rl-Rk), (i, j, Rj-Ri)),
  ]
  append!(model.particle_hole_interactions, Γs)
  append!(model.particle_particle_interactions, Δs)
end
