export HFBComputer
export make_hfbamplitude,
       make_hfbfield,
       compute_hfbfield!,
       make_hamiltonian,
       make_greencollector,
       make_hint,
       randomize!,
       isvalid

export clearΓ!,
       clearΔ!

export make_hoppingmatrix,
       make_Gammamatrix,
       make_Deltamatrix

#using LinearAlgebra
import DataStructures
using HartreeFockBogoliubov

"""
`CollectRow` is holds info on how to compute ρ or t.
Its elements are:
1. Is diagonal? (only for rho)
2. row orbital
3. col orbital
4. displacement r(col) - r(row)
"""
const CollectRow = Tuple{Bool, Int, Int, Vector{Float64}}


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
const DeployRow = Tuple{Bool, Int, Int, Vector{Float64}, Vector{Tuple{Int, ComplexF64, Bool}}}


"""
`HFBConmputer` is a type holding the ρ, t and Γ, Δ of a Hartree-Fock-Bogoliubov
Hamiltonian.
"""
mutable struct HFBComputer{O}
    unitcell ::UnitCell{O}
    hoppings_diagonal ::Vector{Spec.HoppingDiagonal}
    hoppings_offdiagonal ::Vector{Spec.HoppingOffdiagonal}
    temperature ::Float64

    fermi ::Function
    ρ_registry ::Vector{CollectRow}
    t_registry ::Vector{CollectRow}
    Γ_registry ::Vector{DeployRow}
    Δ_registry ::Vector{DeployRow}
end


"""
"""
function HFBComputer(unitcell ::UnitCell{O},
                     hoppings_diagonal ::AbstractVector{Spec.HoppingDiagonal},
                     hoppings_offdiagonal ::AbstractVector{Spec.HoppingOffdiagonal},
                     temperature ::Real) where {O}
    fermi = fermidirac(temperature)
    return HFBComputer{O}(unitcell, hoppings_diagonal, hoppings_offdiagonal, temperature, fermi, [], [], [], [])
end


"""
"""
function make_particleholeregistry(ham::HFBHamiltonian)
    function getdistance(i ::Int, j ::Int, Rij::Vector{Int}) ::Vector{Float64}
        ri, rj = getorbitalcoord(ham.unitcell, i), getorbitalcoord(ham.unitcell, j)
        rj = rj + Rij
        ri, rj = fract2carte(ham.unitcell, ri), fract2carte(ham.unitcell, rj)
        return rj - ri
    end

    collect_reg = Dict()
    deploy_reg = Dict()

    let
        R = zeros(Int, dimension(ham.unitcell))
        r = zeros(Float64, dimension(ham.unitcell))
        for (i, orb) in enumerate(ham.unitcell.orbitals)
            collect_reg[i, i, zero(R)] = (length(collect_reg)+1, (true, i, i, zero(r)))
        end
        for (i, orb) in enumerate(ham.unitcell.orbitals)
            deploy_reg[i, i, zero(R)] = (length(deploy_reg)+1, (true, i, i, zero(r), []))
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
function make_particleparticleregistry(ham::HFBHamiltonian)
    function getdistance(i ::Int, j ::Int, Rij::Vector{Int}) ::Vector{Float64}
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


function HFBComputer(ham ::HFBHamiltonian{O},
                     temperature ::Real;
                     ttol ::Real=eps(Float64),
                     etol ::Real=sqrt(eps(Float64))) where {O}
    dim = dimension(ham.unitcell)

    unitcell = ham.unitcell
    hoppings_diagonal = ham.hoppings_diagonal
    hoppings_offdiagonal = ham.hoppings_offdiagonal
    fermi = fermidirac(temperature)

    (ρ_registry, Γ_registry) = make_particleholeregistry(ham)
    (t_registry, Δ_registry) = make_particleparticleregistry(ham)

    return HFBComputer{O}(unitcell,
                          hoppings_diagonal,
                          hoppings_offdiagonal,
                          temperature, fermi,
                          ρ_registry, t_registry,
                          Γ_registry, Δ_registry)
end


function clearΓ!(computer::HFBComputer)
    for (i, item) in enumerate(computer.Γ_registry)
        computer.Γ_registry[i] = (item[1:end-1]..., [])
    end
end


function clearΔ!(computer::HFBComputer)
    for (i, item) in enumerate(computer.Δ_registry)
        computer.Δ_registry[i] = (item[1:end-1]..., [])
    end
end


"""
func : (idx, i, j, r) -> val
"""
function make_hfbamplitude(computer ::HFBComputer,
                          funcρ ::Function, funct ::Function) ::HFBAmplitude
    ρs = zeros(ComplexF64, length(computer.ρ_registry))
    ts = zeros(ComplexF64, length(computer.t_registry))

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
    return HFBAmplitude(ρs, ts)
end


"""
func : (idx, i, j, r) -> 0
"""
function make_hfbamplitude(computer ::HFBComputer) ::HFBAmplitude
    ρs = zeros(ComplexF64, length(computer.ρ_registry))
    ts = zeros(ComplexF64, length(computer.t_registry))
    return HFBAmplitude(ρs, ts)
end


"""
Compute Γ and Δ from ρ and t.
"""
function make_hfbfield(computer ::HFBComputer) ::HFBField
    Γs = zeros(ComplexF64, length(computer.Γ_registry))
    Δs = zeros(ComplexF64, length(computer.Δ_registry))
    return HFBField(Γs, Δs)
end


"""
Compute Γ and Δ from ρ and t.
"""
function make_hfbfield(computer ::HFBComputer,
                       hfbamplitude ::HFBAmplitude) ::HFBField
    ρs = hfbamplitude.ρ
    ts = hfbamplitude.t
    Γs = zeros(ComplexF64, length(computer.Γ_registry))
    Δs = zeros(ComplexF64, length(computer.Δ_registry))
    for (tgtidx, (_, i, j, r, srcs)) in enumerate(computer.Γ_registry)
        value = zero(ComplexF64)
        for (srcidx, amplitude, star) in srcs
            value += amplitude * (star ? conj(ρs[srcidx]) : ρs[srcidx])
        end
        Γs[tgtidx] = value
    end

    for (tgtidx, (_, i, j, r, srcs)) in enumerate(computer.Δ_registry)
        value = zero(ComplexF64)
        for (srcidx, amplitude, neg) in srcs
            value += amplitude * (neg ? -ts[srcidx] : ts[srcidx])
        end
        Δs[tgtidx] = value
    end
    return HFBField(Γs, Δs)
end


"""
Compute Γ and Δ from ρ and t.
"""
function compute_hfbfield!(hfbfield ::HFBField,
                           computer ::HFBComputer,
                           hfbamplitude ::HFBAmplitude) ::HFBField
    fill!(hfbfield.Γ, 0)
    fill!(hfbfield.Δ, 0)
    ρs = hfbamplitude.ρ
    ts = hfbamplitude.t
    for (tgtidx, (_, i, j, r, srcs)) in enumerate(computer.Γ_registry)
        value = zero(ComplexF64)
        for (srcidx, amplitude, star) in srcs
            value += amplitude * (star ? conj(ρs[srcidx]) : ρs[srcidx])
        end
        hfbfield.Γ[tgtidx] = value
    end

    for (tgtidx, (_, i, j, r, srcs)) in enumerate(computer.Δ_registry)
        value = zero(ComplexF64)
        for (srcidx, amplitude, neg) in srcs
            value += amplitude * (neg ? -ts[srcidx] : ts[srcidx])
        end
        hfbfield.Δ[tgtidx] = value
    end
    hfbfield
end



"""
Return a generator of hopping matrix (which is a function of momentum)
"""
function make_hoppingmatrix(computer::HFBComputer) ::Function
    norb = numorbital(computer.unitcell)
    hk = Generator.hopping_inplace(computer.unitcell, computer.hoppings_diagonal, computer.hoppings_offdiagonal)

    function ret(momentum ::AbstractVector{Float64}) ::Matrix{ComplexF64}
        out = zeros(ComplexF64, (norb, norb))
        hk(momentum, out)
        return out
    end
end


"""
Return a generator of Γ matrix (which is a function of momentum)
"""
function make_Gammamatrix(computer::HFBComputer,
                          theΓs ::AbstractVector{<:Number}) ::Function
    norb = numorbital(computer.unitcell)
    Γ_registry = deepcopy(computer.Γ_registry)
    Γs = Vector{ComplexF64}(theΓs)
    @assert(length(Γ_registry) == length(Γs))
    for ((isdiag, i, j, r, _), Γ) in zip(Γ_registry, Γs)
        if isdiag
            @assert(i==j && all(x -> isapprox(x, 0), r))
            @assert(isapprox(imag(Γ), 0))
        else
            @assert(!(i==j && all(x -> isapprox(x, 0), r)))
        end
    end
    function ret(momentum ::AbstractVector{Float64}) ::Matrix{ComplexF64}
        out = zeros(ComplexF64, (norb, norb))
        for ((isdiag, i, j, r, _), Γ) in zip(Γ_registry, Γs)
            if isdiag
                out[i,i] += real(Γ)
            else
                phase = cis(dot(momentum, r))
                Γp = Γ * phase
                out[i,j] += Γp
                out[j,i] += conj(Γp)
            end
        end
        return out
    end
end

function make_Gammamatrix(computer ::HFBComputer,
                          hfbfield ::HFBField) ::Function
    return make_Gammamatrix(computer, hfbfield.Γ)
end


"""
Return a generator of Δ matrix (which is a function of momentum)
"""
function make_Deltamatrix(computer::HFBComputer,
                          theΔs ::AbstractVector{<:Number}) ::Function
    norb = numorbital(computer.unitcell)
    Δ_registry = deepcopy(computer.Δ_registry)
    Δs = Vector{ComplexF64}(theΔs)
    @assert(length(Δ_registry) == length(Δs))
    for ((isdiag, i, j, r, _), Δ) in zip(Δ_registry, Δs)
        @assert( !isdiag )
        @assert( !(i==j && all(x -> isapprox(x, 0), r)) )
    end
    function ret(momentum ::AbstractVector{Float64}) ::Matrix{ComplexF64}
        out = zeros(ComplexF64, (norb, norb))
        for ((isdiag, i, j, r, _), Δ) in zip(Δ_registry, Δs)
            phase = cis(dot(momentum, r))
            out[i,j] += Δ * phase
            out[j,i] -= Δ * conj(phase)
        end
        return out
    end
end

function make_Deltamatrix(computer ::HFBComputer,
                          hfbfield ::HFBField) ::Function
    return make_Deltamatrix(computer, hfbfield.Δ)
end


"""
"""
function make_hamiltonian(computer ::HFBComputer,
                          theΓs ::AbstractVector{<:Number},
                          theΔs ::AbstractVector{<:Number}) ::Function
    norb = numorbital(computer.unitcell)
    hk = Generator.hopping_inplace(computer.unitcell, computer.hoppings_diagonal, computer.hoppings_offdiagonal)
    Γ_registry = deepcopy(computer.Γ_registry)
    Δ_registry = deepcopy(computer.Δ_registry)
    Γs = Vector{ComplexF64}(theΓs)
    Δs = Vector{ComplexF64}(theΔs)

    @assert(length(Γ_registry) == length(Γs))
    @assert(length(Δ_registry) == length(Δs))

    for ((isdiag, i, j, r, _), Γ) in zip(Γ_registry, Γs)
        if isdiag
            @assert(i==j && all(x -> isapprox(x, 0), r))
            @assert(isapprox(imag(Γ), 0))
        else
            @assert(!(i==j && all(x -> isapprox(x, 0), r)))
        end
    end
    for ((isdiag, i, j, r, _), Δ) in zip(Δ_registry, Δs)
        @assert( !isdiag )
        @assert( !(i==j && all(x -> isapprox(x, 0), r)) )
    end

    function ret(momentum ::AbstractVector{Float64}) ::Matrix{ComplexF64}
        #out = zeros(ComplexF64, (norb, 2, norb, 2))
        #h11 = view(out, :, 1, :, 1)
        #h12 = view(out, :, 1, :, 2)
        #h21 = view(out, :, 2, :, 1)
        #h22 = view(out, :, 2, :, 2)
        out = zeros(ComplexF64, (2, norb, 2, norb))
        h11 = view(out, 1, :, 1, :)
        h12 = view(out, 1, :, 2, :)
        h21 = view(out, 2, :, 1, :)
        h22 = view(out, 2, :, 2, :)

        # 1/3. non-interacting kinetic part
        hk( momentum, h11)
        hk(-momentum, h22)

        # 2/3. Gamma
        for ((isdiag, i, j, r, _), Γ) in zip(Γ_registry, Γs)
            if isdiag
                h11[i,i] += real(Γ)
                h22[i,i] += real(Γ) # Note: The sign will be flipped after the loop
            else
                phase = cis(dot(momentum, r))
                h11[i,j] += Γ * phase
                h22[i,j] += Γ * conj(phase)
                h11[j,i] += conj(Γ * phase)
                h22[j,i] += conj(Γ) * phase
            end
        end
        h22[:,:] = -transpose(h22)

        # 3/3. Delta
        for ((isdiag, i, j, r, _), Δ) in zip(Δ_registry, Δs)
            phase = cis(dot(momentum, r))
            val1 = Δ * phase
            val2 = Δ * conj(phase)
            h12[i,j] += val1
            h12[j,i] -= val2
            h21[j,i] += conj(val1)
            h21[i,j] -= conj(val2)
        end
        return reshape(out, (2*norb, 2*norb))
    end
end

function make_hamiltonian(computer ::HFBComputer,
                          hfbfield ::HFBField) ::Function
    return make_hamiltonian(computer, hfbfield.Γ, hfbfield.Δ)
end

function make_hamiltonian(computer ::HFBComputer,
                          hfbamplitude ::HFBAmplitude) ::Function
    hfbfield = make_hfbfield(computer, hfbamplitude)
    return make_hamiltonian(computer, hfbfield.Γ, hfbfield.Δ)
end


"""
    make_greencollector

Returns a function which has the following signature
```
collector(k, eigenvalues, eigenvectors, ρout, tout)
```
"""
function make_greencollector(computer::HFBComputer{O}) ::Function where {O}
    fermi = computer.fermi
    norb = numorbital(computer.unitcell)
    ρ_registry = deepcopy(computer.ρ_registry)
    t_registry = deepcopy(computer.t_registry)

    for (isdiag, i, j, r) in ρ_registry
        @assert(isdiag == (i==j && all(x -> isapprox(x, 0), r)))
    end
    for (isdiag, i, j, r) in t_registry
        @assert(!isdiag && !(i==j && all(x -> isapprox(x, 0), r)))
    end

    function ret(k::AbstractVector{<:Real},
                 eigenvalues ::AbstractVector{<:Real},
                 eigenvectors ::AbstractMatrix{<:Number},
                 hfbamplitude ::HFBAmplitude) ::Nothing

        @assert(length(eigenvalues) == 2*norb)
        @assert(size(eigenvectors) == (2*norb, 2*norb))

        f = [fermi(e) for e in eigenvalues]
        ψ = reshape(eigenvectors, (2, norb, 2*norb))
        u = view(ψ, 1, :, :)
        v = view(ψ, 2, :, :)

        uf = u * Diagonal(f)
        ρmat = uf * (u')
        tmat = uf * (v')

        for (idx, (isdiag, i, j, r)) in enumerate(ρ_registry)
            if isdiag
                hfbamplitude.ρ[idx] += real( ρmat[i,i] )
            else
                hfbamplitude.ρ[idx] += ρmat[i, j] * cis(-dot(k, r))
            end
        end
        for (idx, (isdiag, i, j, r)) in enumerate(t_registry)
            hfbamplitude.t[idx] += tmat[i, j] * cis(-dot(k, r))
        end
        nothing
    end
end


"""
    Check if hint contains ρ
"""
function make_hfbamplitude(computer::HFBComputer{O}, hint::HFBAmplitudeHint{O}) ::HFBAmplitude where {O}
    hfbamplitude = HFB.make_hfbamplitude(computer)
    unitcell = computer.unitcell

    for (ρidx, (i, j, r)) in enumerate(computer.ρ_registry)
        iname, icoord = getorbital(unitcell, i)
        jname, jcoord = getorbital(unitcell, j)
        Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
        Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)
        R = Rj - Ri
        if haskey(hint.ρ, (iname, jname, R))
            hfbamplitude.ρ[ρidx] = hint.ρ[iname, jname, R]
        elseif haskey(hint.ρ, (jname, iname, -R))
            hfbamplitude.ρ[ρidx] = conj(hint.ρ[jname, iname, -R])
        end
    end

    for (tidx, (i, j, r)) in enumerate(computer.t_registry)
        iname, icoord = getorbital(unitcell, i)
        jname, jcoord = getorbital(unitcell, j)
        Ri = whichunitcell(unitcell, iname, fract2carte(unitcell, icoord))
        Rj = whichunitcell(unitcell, jname, fract2carte(unitcell, icoord) + r)
        R = Rj - Ri
        if haskey(hint.t, (iname, jname, R))
            hfbamplitude.t[tidx] = hint.t[iname, jname, R]
        elseif haskey(hint.t, (jname, iname, -R))
            hfbamplitude.t[tidx] = conj(hint.t[jname, iname, -R])
        end
    end
    return hfbamplitude
end


"""
"""
function make_hint(computer::HFBComputer{O}, sol::HFBAmplitude) ::HFBAmplitudeHint{O} where {O}
    ρ = DataStructures.OrderedDict{Tuple{O, O, Vector{Int}}, ComplexF64}()
    t = DataStructures.OrderedDict{Tuple{O, O, Vector{Int}}, ComplexF64}()
    uc = computer.unitcell

    for (ρidx, (_, i, j, rij)) in enumerate(computer.ρ_registry)
        iname, icoord = getorbital(uc, i)
        jname, jcoord = getorbital(uc, j)
        Ri = whichunitcell(uc, iname, fract2carte(uc, icoord))
        Rj = whichunitcell(uc, jname, fract2carte(uc, icoord) + rij)
        @assert(!haskey(ρ, (iname, jname, Rj-Ri)))
        ρ[iname,jname,Rj-Ri] = sol.ρs[ρidx]
    end

    for (tidx, (_, i, j, rij)) in enumerate(computer.t_registry)
        iname, icoord = getorbital(uc, i)
        jname, jcoord = getorbital(uc, j)
        Ri = whichunitcell(uc, iname, fract2carte(uc, icoord))
        Rj = whichunitcell(uc, jname, fract2carte(uc, icoord) + rij)
        @assert(!haskey(t, (iname, jname, Rj-Ri)))
        t[iname,jname,Rj-Ri] = sol.ts[tidx]
    end

    return HFBAmplitudeHint{O}(ρ, t)
end


"""
Randomize a hfbamplitude
"""
function randomize!(computer::HFBComputer, sol ::HFBAmplitude)
    sol.ρ[:] = rand(Float64, length(sol.ρ))
    sol.t[:] = rand(ComplexF64, length(sol.t))
end


"""
"""
function isvalid(computer ::HFBComputer,
                 solution ::HFBAmplitude;
                 tolerance ::Real=sqrt(eps(Float64))) ::Bool
    if (length(computer.ρ_registry) != length(solution.ρ) ||
        length(computer.t_registry) != length(solution.t))
        return false
    end
    for (idx, (isdiag, i, j, r)) in enumerate(computer.ρ_registry)
        if isdiag && ! isapprox(imag(solution.ρ[idx]), 0; atol=tolerance)
            return false
        end
    end
    return true
end


"""
"""
function isvalid(computer ::HFBComputer,
                 solution ::HFBField;
                 tolerance ::Real=sqrt(eps(Float64))) ::Bool
    if (length(computer.Γ_registry) != length(solution.Γ) ||
        length(computer.Δ_registry) != length(solution.Δ) )
        return false
    end
    for (idx, (isdiag, i, j, r, _)) in enumerate(computer.Γ_registry)
        if isdiag && ! isapprox(imag(solution.Γ[idx]), 0; atol=tolerance)
            return false
        end
    end
    return true
end
