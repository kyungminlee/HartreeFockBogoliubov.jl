#export addinteraction!

export HFBHamiltonian
export addinteraction!


struct HoppingMeanField{R<:Real}
  amplitude ::R
  target ::Tuple{Int64, Int64, Vector{Int64}}
  source ::Tuple{Int64, Int64, Vector{Int64}}
  star ::Bool

  function HoppingMeanField{R}(
          amplitude ::R,
          target::Tuple{Int64, Int64, Vector{Int64}},
          source::Tuple{Int64, Int64, Vector{Int64}},
          star ::Bool
        ) where {R<:Real}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    @assert(i <= j)
    @assert(k <= l)
    new{R}(amplitude, target, source, star)
  end

  function HoppingMeanField(
          amplitude ::C,
          target::Tuple{Int64, Int64, Vector{Int64}},
          source::Tuple{Int64, Int64, Vector{Int64}}
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
    new{R}(amplitude, (i,j,Rij), (k,l,Rkl), star)
  end
end


struct PairingMeanField{R<:Real}
  amplitude ::R
  target ::Tuple{Int64, Int64, Vector{Int64}}
  source ::Tuple{Int64, Int64, Vector{Int64}}
  negate ::Bool

  function PairingMeanField{R}(amplitude ::R,
                               target::Tuple{Int64, Int64, Vector{Int64}},
                               source::Tuple{Int64, Int64, Vector{Int64}},
                               negate ::Bool
                               ) where {R<:Real}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    @assert(i <= j)
    @assert(k <= l)
    new{R}(amplitude, target, source, negate)
  end

  function PairingMeanField{R}(
          amplitude ::R,
          target::Tuple{Int64, Int64, Vector{Int64}},
          source::Tuple{Int64, Int64, Vector{Int64}}
        ) where {R<:Real}
    (i, j, Rij) = target
    (k, l, Rkl) = source
    negate = false
    if i > j
      (i,j) = (j,i)
      Rij = -Rij
      negate = !negate
    end
    if k > l
      (k,l) = (l,k)
      Rkl = -Rkl
      negate = !negate
    end
    new{R}(amplitude, (i,j,Rij), (k,l,Rkl), negate)
  end
end


mutable struct HFBHamiltonian{O}
  unitcell ::UnitCell{O}
  hoppings ::Vector{SpecHopping}
  particle_hole_interactions ::Vector{HoppingMeanField}
  particle_particle_interactions ::Vector{PairingMeanField}
end


function HFBHamiltonian(
        unitcell ::UnitCell{O},
        hoppings ::AbstractVector{SpecHopping},
        particle_hole_interactions ::AbstractVector{HoppingMeanField},
        particle_particle_interactions ::AbstractVector{PairingMeanField}
      ) where {O}
  return HFBHamiltonian{O}(unitcell,
                           hoppings,
                           particle_hole_interactions,
                           particle_particle_interactions)
end


function HFBHamiltonian(full_hamiltonian ::SpecHamiltonian{O}) where {O}
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
                         specint ::SpecInteractionDiagonal{R}) where {O,R<:Real}
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
