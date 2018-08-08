
struct MeasureParticleHoleDiagonal
    i ::Int
    j ::Int
    Rij ::Vector{Int}
    rij ::Vector{Float64}

    function MeasureParticleHoleDiagonal(
        unitcell::UnitCell,
        i::Integer)
        ndim = dimension(unitcell)
        return new(i, i, zeros(Int, ndim), zeros(Float64, ndim))
    end
end

struct MeasureParticleHoleOffdiagonal
    i ::Int
    j ::Int
    Rij ::Vector{Int}
    rij ::Vector{Float64}

    function MeasureParticleHoleOffdiagonal(
        unitcell ::UnitCell,
        i::Integer,
        j::Integer,
        Rij::AbstractVector{I}) where {I<:Integer}
        @assert( i != j || any(x -> x != 0, Rij) )

        ri, rj = getorbitalcoord(unitcell, i), getorbitalcoord(unitcell, j)
        rj = rj + Rij
        ri, rj = fract2carte(unitcell, ri), fract2carte(unitcell, rj)
        rij = rj - ri

        if i < j
            return new(i,j,Rij, rij)
        else
            return new(j,i,-Rij, -rij)
        end
    end
end

struct MeasureParticleParticle
    i ::Int
    j ::Int
    Rij ::Vector{Int}
    rij ::Vector{Float64}

    function MeasureParticleParticle(
                unitcell::UnitCell,
                i::Integer,
                j::Integer,
                Rij::AbstractVector{I},
            ) where {I<:Integer}
        @assert( i != j || any(x -> x != 0, Rij) )

        ri, rj = getorbitalcoord(unitcell, i), getorbitalcoord(unitcell, j)
        rj = rj + Rij
        ri, rj = fract2carte(unitcell, ri), fract2carte(unitcell, rj)
        rij = rj - ri

        if i < j
            return new(i,j, Rij,  rij)
        else
            return new(j,i,-Rij, -rij)
        end
    end
end

isdiagonal(m ::MeasureParticleHoleDiagonal) = true
isdiagonal(m ::MeasureParticleHoleOffdiagonal) = false
isdiagonal(m ::MeasureParticleParticle) = false

const MeasureParticleHole = Union{MeasureParticleHoleDiagonal, MeasureParticleHoleOffdiagonal}

struct MeasurementList
    particlehole ::Vector{MeasureParticleHole}
    particleparticle ::Vector{MeasureParticleParticle}

    particlehole_index ::Dict{Tuple{Int, Int, Vector{Int}}, Int}
    particleparticle_index ::Dict{Tuple{Int, Int, Vector{Int}}, Int}
end

function addmeasurement!(ml::MeasurementList, ph::MeasureParticleHole)
    if haskey(ml.particlehole_index, (ph.i, ph.j, ph.Rij))
        return # do nothing
    else
        push!(ml.particlehole, ph)
        ml.particlehole_index[ph.i, ph.j, ph.Rij] = length(ml.particlehole)
    end
end


function addmeasurement!(ml::MeasurementList, pp::MeasureParticleParticle)
    if haskey(ml.particleparticle_index, (pp.i, pp.j, pp.Rij))
        return # do nothing
    else
        push!(ml.particlehole, pp)
        ml.particleparticle_index[pp.i, pp.j, pp.Rij] = length(ml.particleparticle)
    end
end



struct MeasurementValue
    particlehole ::Vector{Complex{Float64}}
    particleparticle ::Vector{Complex{Float64}}

    function MeasurementValue(ml::MeasurementList)
        ph = zeros(Complex{Float64}, length(ml.particlehole))
        pp = zeros(Complex{Float64}, length(ml.particleparticle))
        return new(ph, pp)
    end
end

function clear(mv::MeasurementValue)
    mv.particlehole[:] = 0.0
    mv.particleparticle[:] = 0.0
    mv
end

function measure(uc::UnitCell{O},
    ph::MeasureParticleHole) where {O}

end

function measure(uc::UnitCell{O},
    ml::MeasurementList, mv::MeasurementValue,
    ρ::Function, t::Function) where {O}
    #for (idx, ph) in enumerate(ml.particlehole)
    #  mv.particlehole[idx] += ρ(i,j) * cis(rij)
end
