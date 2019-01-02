export HFBSolver
export next_hfbamplitude,
       next_hfbamplitude!,
       #next_hfbsolution,
       #next_hfbfield,
       next_hfbamplitude_threaded,
       next_hfbamplitude_threaded!,
       loop,
       loop_threaded,
       #loop_python,
       simpleupdate,
       isvalid

import LinearAlgebra

function eigh(mat ::Matrix{R}) where {R<:AbstractFloat}
    return (LinearAlgebra.eigen(LinearAlgebra.Symmetric(mat))...,)
end

function eigh(mat ::Matrix{C}) where {C<:Complex}
    return (LinearAlgebra.eigen(LinearAlgebra.Hermitian(mat))...,)
end


"""
"""
mutable struct HFBSolver{O}
    # Originals
    hamiltonian ::FullHamiltonian{O}
    size ::Vector{Int}
    temperature ::Float64

    # Derivatives
    momentumgrid ::Array{Vector{Float64}}
    hfbhamiltonian ::HFBHamiltonian{O}
    hfbcomputer ::HFBComputer{O}
    greencollector ::Function

    # Optional
    eigensolver ::Function
end


"""
"""
function HFBSolver(fullhamiltonian::FullHamiltonian{O},
                   size ::AbstractVector{<:Integer},
                   temperature ::Real;
                   eigensolver::Function=eigh,
                   tol ::Real=eps(Float64)) ::HFBSolver{O} where {O}
    dim = dimension(fullhamiltonian.unitcell)
    if length(size) != dim
        throw(ArgumentError("dimension and size do not match: size=$(size)"))
    elseif !all((x) -> x > 0, size)
        throw(ArgumentError("size should be positive: size=$(size)"))
    elseif temperature < 0
        throw(ArgumentError("temperature should be non-negative"))
    end
    momentumgrid = Lattice.momentumgrid(fullhamiltonian.unitcell, size)
    hfbhamiltonian = HFB.HFBHamiltonian(fullhamiltonian)
    hfbcomputer = HFB.HFBComputer(hfbhamiltonian, temperature)
    greencollector = HFB.make_greencollector(hfbcomputer)
    return HFBSolver{O}(fullhamiltonian, size, temperature,
                        momentumgrid,
                        hfbhamiltonian, hfbcomputer, greencollector,
                        eigensolver)
end


"""
"""
function next_hfbamplitude!(next_hfbamplitude ::HFBAmplitude,
                            solver ::HFBSolver,
                            hfbfield ::HFBField)
    fill!(next_hfbamplitude.ρ, 0)
    fill!(next_hfbamplitude.t, 0)
    ham = make_hamiltonian(solver.hfbcomputer, hfbfield)
    for momentum in solver.momentumgrid
        hk = ham(momentum)
        (eivals, eivecs) = solver.eigensolver(hk)
        #(eivals, eivecs) = (LinearAlgebra.eigen(LinearAlgebra.Hermitian(hk))...,)
        solver.greencollector(momentum, eivals, eivecs, next_hfbamplitude)
    end
    next_hfbamplitude.ρ /= length(solver.momentumgrid)
    next_hfbamplitude.t /= length(solver.momentumgrid)
    next_hfbamplitude
end


"""
"""
function next_hfbamplitude(solver ::HFBSolver,
                           hfbfield ::HFBField) ::HFBAmplitude
    next_hfbamplitude = make_hfbamplitude(solver.hfbcomputer)
    ham = make_hamiltonian(solver.hfbcomputer, hfbfield)
    for momentum in solver.momentumgrid
        hk = ham(momentum)
        #(eivals, eivecs) = (LinearAlgebra.eigen(LinearAlgebra.Hermitian(hk))...,)
        (eivals, eivecs) = solver.eigensolver(hk)
        solver.greencollector(momentum, eivals, eivecs, next_hfbamplitude)
    end
    next_hfbamplitude.ρ /= length(solver.momentumgrid)
    next_hfbamplitude.t /= length(solver.momentumgrid)
    return next_hfbamplitude
end


"""
"""
function next_hfbamplitude(solver ::HFBSolver,
                           hfbamplitude ::HFBAmplitude) ::HFBAmplitude
    hfbfield = make_hfbfield(solver.hfbcomputer, hfbamplitude)
    return next_hfbamplitude(solver, hfbfield)
end


"""
"""
function next_hfbfield(solver ::HFBSolver, hfbfield ::HFBField) ::HFBAmplitude
    nsf = make_hfbamplitude(solver.hfbcomputer)
    ham = make_hamiltonian(solver.hfbcomputer, hfbfield)
    for momentum in solver.momentumgrid
        hk = ham(momentum)
        #(eivals, eivecs) = (LinearAlgebra.eigen(LinearAlgebra.Hermitian(hk))...,)
        (eivals, eivecs) = solver.eigensolver(hk)
        solver.greencollector(momentum, eivals, eivecs, nsf)
    end
    nsf.ρ /= length(solver.momentumgrid)
    nsf.t /= length(solver.momentumgrid)
    return make_hfbfield(solver, nsf)
end


#=
"""
"""
function next_hfbsolution(solver ::HFBSolver, hfbfield ::HFBField) ::HFBAmplitude
    nsf = make_hfbamplitude(solver.hfbcomputer)
    ham = make_hamiltonian(solver.hfbcomputer, hfbfield)
    for momentum in solver.momentumgrid
        hk = ham(momentum)
        #(eivals, eivecs) = (LinearAlgebra.eigen(LinearAlgebra.Hermitian(hk))...,)
        (eivals, eivecs) = solver.eigensolver(hk)

        solver.greencollector(momentum, eivals, eivecs, nsf)
    end
    nsf.ρ /= length(solver.momentumgrid)
    nsf.t /= length(solver.momentumgrid)
    return (nsf, make_hfbfield(solver, nsf))
end
=#


#=
"""
"""
function getnextsolutionpython(solver ::HFBSolver{T}, sol ::HFBAmplitude) ::HFBAmplitude where {T}
    newsol = make_hfbfield(solver.hfbcomputer)
    ham = makehamiltonian(solver.hfbcomputer, sol.Γ, sol.Δ)
    for momentum in solver.momentumgrid
        hk = ham(momentum)
        (eivals, eivecs) = npl[:eigh](hk)
        solver.greencollector(momentum, eivals, eivecs, newsol.ρ, newsol.t)
    end
    newsol.ρ /= length(solver.momentumgrid)
    newsol.t /= length(solver.momentumgrid)
    newsol.Γ[:], newsol.Δ[:] = HFB.computehfbfields(solver.hfbcomputer, newsol.ρ, newsol.t)
    return newsol
end
=#

"""
"""
function next_hfbamplitude_threaded(solver ::HFBSolver,
                                   hfbfield ::HFBField) ::HFBAmplitude
    nsf = make_hfbamplitude(solver.hfbcomputer)
    ham = make_hamiltonian(solver.hfbcomputer, hfbfield)

    #threadedρ = [zero(nsf.ρ) for tid in 1:Threads.nthreads()]
    #threadedt = [zero(nsf.t) for tid in 1:Threads.nthreads()]
    threaded_nsf = [make_hfbamplitude(solver) for tid in 1:Threads.nthreads()]

    num_momentum = length(solver.momentumgrid)
    Threads.@threads for idx_momentum in 1:num_momentum
        momentum = solver.momentumgrid[idx_momentum]
        tid = Threads.threadid()
        hk = ham(momentum)
        #(eivals, eivecs) = (LinearAlgebra.eigen(LinearAlgebra.Hermitian(hk))...,)
        (eivals, eivecs) = solver.eigensolver(hk)
        solver.greencollector(momentum, eivals, eivecs, threaded_nsf[tid])
    end

    for tnsf in threaded_nsf
        nsf.ρ += tnsf.ρ
        nsf.t += tnsf.t
    end

    nsf.ρ /= length(solver.momentumgrid)
    nsf.t /= length(solver.momentumgrid)
    return nsf
end



"""
"""
function next_hfbamplitude_threaded!(next_hfbamplitude ::HFBAmplitude,
                                     solver ::HFBSolver,
                                     hfbfield ::HFBField)
    fill!(next_hfbamplitude.ρ, 0)
    fill!(next_hfbamplitude.t, 0)
    ham = make_hamiltonian(solver.hfbcomputer, hfbfield)

    threaded_nsf = [make_hfbamplitude(solver) for tid in 1:Threads.nthreads()]

    num_momentum = length(solver.momentumgrid)
    Threads.@threads for idx_momentum in 1:num_momentum
        momentum = solver.momentumgrid[idx_momentum]
        tid = Threads.threadid()
        hk = ham(momentum)
        #(eivals, eivecs) = (LinearAlgebra.eigen(LinearAlgebra.Hermitian(hk))...,)
        (eivals, eivecs) = solver.eigensolver(hk)
        solver.greencollector(momentum, eivals, eivecs, threaded_nsf[tid])
    end

    for tnsf in threaded_nsf
        next_hfbamplitude.ρ += tnsf.ρ
        next_hfbamplitude.t += tnsf.t
    end

    next_hfbamplitude.ρ /= length(solver.momentumgrid)
    next_hfbamplitude.t /= length(solver.momentumgrid)
    return next_hfbamplitude
end


function simpleupdate(sol::HFBField, newsol ::HFBField) ::HFBField
    # Tested anyway by the assignment
    sol.Γ[:] = newsol.Γ[:]
    sol.Δ[:] = newsol.Δ[:]
    sol
end



_noop(x...) = nothing


"""
    loop

Perform selfconsistency loop a number of times with the given precondition and given update functions.

# Arguments
* `solver ::HFBSolver{T}`
* `sol::HFBAmplitude`
* `run::Integer`

# Optional Arguments
* `update::Function=simpleupdate`
* `precondition::Function=identity`: on target field
* `callback::Function=_noop`: Function called after every update as
`callback(i, run)`
"""
function loop(solver ::HFBSolver,
              hfbfield ::HFBField,
              run::Integer;
              update::Function=simpleupdate,
              precondition::Function=identity,
              callback::Function=_noop)
    hfbamplitude = make_hfbamplitude(solver)
    new_hfbfield = make_hfbfield(solver)
    for i in 1:run
        precondition(hfbfield)
        next_hfbamplitude!(hfbamplitude, solver, hfbfield)
        compute_hfbfield!(new_hfbfield, solver, hfbamplitude)
        update(hfbfield, new_hfbfield)
        callback(i, run)
    end
    (hfbamplitude, hfbfield)
end


"""
    loop

Perform selfconsistency loop a number of times with the given precondition and given update functions.

# Arguments
* `solver ::HFBSolver{T}`
* `sf::HFBAmplitude`
* `run::Integer`

# Optional Arguments
* `update::Function=simpleupdate`
* `precondition::Function=identity`
* `callback::Function=_noop`: Function called after every update as
`callback(i, run)`
"""
function loop_threaded(solver ::HFBSolver,
                       hfbamplitude ::HFBAmplitude,
                       run ::Integer;
                       update::Function=simpleupdate,
                       precondition::Function=identity,
                       callback::Function=_noop)
    #hfbfield = make_hfbfield(solver, hfbamplitude)
    hfbamplitude = nothing
    for i in 1:run
        precondition(hfbfield)
        next_hfbamplitude_threaded!(hfbamplitude, solver, hfbfield)
        update(hfbfield, make_hfbfield(solver, hfbamplitude))
        callback(i, run)
    end
    (hfbamplitude, hfbfield)
end


#=
"""
    loop

Perform selfconsistency loop a number of times with the given precondition and given update functions
using Python's `numpy.linalg.eigh` (which hopefully is using MKL library).

# Arguments
* `solver ::HFBSolver{T}`
* `sf::HFBAmplitude`
* `run::Integer`

# Optional Arguments
* `update::Function=simpleupdate`
* `precondition::Function=identity`
* `callback::Function=_noop`: Function called after every update as
`callback(i, run)`
"""
function looppython(solver ::HFBSolver{O},
                    sf::HFBAmplitude,
                    run::Integer;
                    update::Function=simpleupdate,
                    precondition::Function=identity,
                    callback::Function=_noop,
                    ) ::HFBAmplitude where {O}
    sol = copy(sol)
    for i in 1:run
        precondition(sol)
        newsol = getnextsolutionpython(solver, sol)
        update(sol, newsol)
        callback(i, run)
    end
    sol
end
=#


# ===== Functions from HFBComputer
for fname in [:make_hfbamplitude,
              :make_hfbfield,
              :make_hoppingmatrix, :make_Gammamatrix, :make_Deltamatrix,
              :make_hamiltonian,
              :make_greencollector,
              :make_hint,
              :randomize!,
              :freeze,
              :isvalid
              ]
    @eval begin
        $fname(solver::HFBSolver, args...) = ($fname)(solver.hfbcomputer, args...)
    end
end
