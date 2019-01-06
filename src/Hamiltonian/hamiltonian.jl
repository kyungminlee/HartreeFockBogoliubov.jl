"""
    FullHamiltonian

  # Members
  * `unitcell ::UnitCell`
  * `hoppings ::Vector{Hopping}`
  * `interactions ::Vector{Interaction}`
"""
mutable struct FullHamiltonian{O}
    unitcell ::UnitCell{O}
    hoppings_diagonal ::Vector{HoppingDiagonal}
    hoppings_offdiagonal ::Vector{HoppingOffdiagonal}
    interactions_diagonal ::Vector{InteractionDiagonal}
    interactions_offdiagonal ::Vector{InteractionOffdiagonal}
end


"""
    Hamiltonian

  Create an empty Hamiltonian

  # Arguments
  * `unitcell ::UnitCell`
"""
function FullHamiltonian(unitcell::UnitCell{O}) where {O}
    return FullHamiltonian{O}(unitcell,
                              HoppingDiagonal[],
                              HoppingOffdiagonal[],
                              InteractionDiagonal[],
                              InteractionOffdiagonal[])
end

"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingDiagonal`
"""
function addhopping!(hamiltonian ::FullHamiltonian{O},
                     hopping ::HoppingDiagonal{R};
                     tol::Real=eps(Float64)) where {O, R<:Real}
    @assert(1 <= hopping.i <= numorbital(hamiltonian.unitcell),
            "orbital $(hopping.i) not defined in unit cell")
    if abs(hopping.amplitude) > tol
        push!(hamiltonian.hoppings_diagonal, hopping)
    end
end


"""
    addhopping!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `hopping ::HoppingOffdiagonal`
"""
function addhopping!(hamiltonian ::FullHamiltonian{O},
                     hopping ::HoppingOffdiagonal{C};
                     tol ::Real=eps(Float64)) where {O, C<:Number}
    @assert(1 <= hopping.j <= numorbital(hamiltonian.unitcell),
            "orbital $(hopping.i) not defined in unit cell")
    @assert(1 <= hopping.j <= numorbital(hamiltonian.unitcell),
            "orbital $(hopping.j) not defined in unit cell")
    if abs(hopping.amplitude) > tol
        push!(hamiltonian.hoppings_offdiagonal, hopping)
    end
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionDiagonal`
"""
function addinteraction!(hamiltonian ::FullHamiltonian{O},
                         interaction ::InteractionDiagonal{R};
                         tol ::Real=eps(Float64)) where {O, R<:Real}
    @assert(1 <= interaction.i <= numorbital(hamiltonian.unitcell),
            "orbital $(interaction.i) not defined in unit cell")
    @assert(1 <= interaction.j <= numorbital(hamiltonian.unitcell),
            "orbital $(interaction.j) not defined in unit cell")
    if abs(interaction.amplitude) > tol
        push!(hamiltonian.interactions_diagonal, interaction)
    end
end


"""
    addinteraction!

  # Arguments
  * `hamiltonian ::Hamiltonian`
  * `interaction ::InteractionOffdiagonal`
"""
function addinteraction!(hamiltonian ::FullHamiltonian{O},
                         interaction::InteractionOffdiagonal{C};
                         tol=eps(Float64)) where {O, C<:Number}
    @assert(1 <= interaction.i <= numorbital(hamiltonian.unitcell),
            "orbital $(interaction.i) not defined in unit cell")
    @assert(1 <= interaction.j <= numorbital(hamiltonian.unitcell),
            "orbital $(interaction.j) not defined in unit cell")
    @assert(1 <= interaction.k <= numorbital(hamiltonian.unitcell),
            "orbital $(interaction.k) not defined in unit cell")
    @assert(1 <= interaction.l <= numorbital(hamiltonian.unitcell),
            "orbital $(interaction.l) not defined in unit cell")
    if abs(interaction.amplitude) > tol
        push!(hamiltonian.interactions_offdiagonal, interaction)
    end
end


function simplify(hamiltonian::FullHamiltonian{O}) where {O}
    hopdiadict = Dict{Int, Float64}()
    hopoffdict = Dict{Tuple{Int, Int, Vector{Int}}, ComplexF64}()

    #=
    for hop in hoppings
        if isa(newhop, HoppingDiagonal)
            k = hop.i
            if haskey(hopdiadict, k)
                hopdiadict[k].amplitude += hop.amplitude
            else
                hopdiadict[k] = hop
            end
        else
            k = (hop.i, hop.j, hop.Rj - hop.Ri)
            if haskey(hopoffdict, k)
                hopoffdict[k].amplitude += hop.amplitude
            else
                hopoffdict[k] = hop
            end
        end
    end
    =#

    for hop in hoppings_diagonal
            k = hop.i
            if haskey(hopdiadict, k)
                hopdiadict[k].amplitude += hop.amplitude
            else
                hopdiadict[k] = hop
            end
    end
    for hop in hoppings_offdiagonal
        k = (hop.i, hop.j, hop.Rj - hop.Ri)
        if haskey(hopoffdict, k)
            hopoffdict[k].amplitude += hop.amplitude
        else
            hopoffdict[k] = hop
        end
    end

    # TODO
    error("unimplemented")

    #newhamiltonian = FullHamiltonian{O}(hamiltonian.unitcell)
    #for hop in hamiltonian.hoppings
    #    addhopping!
    #end
end
