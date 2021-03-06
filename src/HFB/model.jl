export HFBHamiltonian
export addinteraction!


"""
"""
struct HoppingMeanField{C<:Number}
    amplitude ::C
    target ::Tuple{Int, Int, Vector{Int}}
    source ::Tuple{Int, Int, Vector{Int}}
    star ::Bool
    function HoppingMeanField{C}(
            amplitude ::C,
            target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
            source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
            star ::Bool,
            ) where {C<:Number}
        (i, j, Rij) = target
        (k, l, Rkl) = source
        @assert(i <= j)
        @assert(k <= l)
        new{C}(amplitude, target, source, star)
    end
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
    HoppingMeanField{C}(amplitude, (i,j,Rij), (k,l,Rkl), star)
end


"""
"""
struct PairingMeanField{C<:Number}
    amplitude ::C
    target ::Tuple{Int, Int, Vector{Int}}
    source ::Tuple{Int, Int, Vector{Int}}
    negate ::Bool
    function PairingMeanField{C}(
                amplitude ::C,
                target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
                source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
                negate ::Bool,
            ) where {C<:Number}
        (i, j, Rij) = target
        (k, l, Rkl) = source
        @assert(i <= j)
        @assert(k <= l)
        new{C}(amplitude, target, source, negate)
  end
end


function PairingMeanField(
            amplitude ::C,
            target::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
            source::Tuple{<:Integer, <:Integer, <:AbstractVector{<:Integer}},
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
    PairingMeanField{C}(amplitude, (i,j,Rij), (k,l,Rkl), negate)
end


"""
"""
mutable struct HFBHamiltonian{O}
    unitcell ::UnitCell{O}
    hoppings_diagonal ::Vector{HoppingDiagonal}
    hoppings_offdiagonal ::Vector{HoppingOffdiagonal}
    particle_hole_interactions ::Vector{HoppingMeanField}
    particle_particle_interactions ::Vector{PairingMeanField}
end


"""
"""
function HFBHamiltonian(
            unitcell ::UnitCell{O},
            hoppings_diagonal ::AbstractVector{HoppingDiagonal},
            hoppings_offdiagonal ::AbstractVector{HoppingOffdiagonal},
            particle_hole_interactions ::AbstractVector{HoppingMeanField},
            particle_particle_interactions ::AbstractVector{PairingMeanField},
        ) where {O}
    return HFBHamiltonian{O}(unitcell,
                             hoppings_diagonal,
                             hoppings_offdiagonal,
                             particle_hole_interactions,
                             particle_particle_interactions)
end


"""
"""
function HFBHamiltonian(fullhamiltonian ::FullHamiltonian{O}) where {O}
    unitcell = fullhamiltonian.unitcell
    hoppings_diagonal = fullhamiltonian.hoppings_diagonal
    hoppings_offdiagonal = fullhamiltonian.hoppings_offdiagonal
    model = HFBHamiltonian{O}(unitcell, hoppings_diagonal, hoppings_offdiagonal, [], [])
    for interaction in fullhamiltonian.interactions_diagonal
        addinteraction!(model, interaction)
    end
    for interaction in fullhamiltonian.interactions_offdiagonal
        addinteraction!(model, interaction)
    end
    return model
end


"""
    addinteraction!

Add diagonal interaction
"""
function addinteraction!(model ::HFBHamiltonian{O},
                         specint ::InteractionDiagonal{R}) where {O,R<:Real}
    v = specint.amplitude
    (i, j) = (specint.i, specint.j)
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
Add offdiagonal interaction
"""
function addinteraction!(model ::HFBHamiltonian{O},
                         specint ::InteractionOffdiagonal{C}) where {O,C<:Number}
    v = specint.amplitude
    (i, j) = (specint.i, specint.j)
    (k, l) = (specint.k, specint.l)
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
