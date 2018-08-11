using JSON
using DataStructures
using ArgParse
using MicroLogging

using ProgressMeter
using HartreeFockBogoliubov
import HartreeFockBogoliubov: Spec, Generator, HFB
using HartreeFockBogoliubov: HFB

function detectresult(outpath::AbstractString)
    for resultfilename in ["result.json",]
        result = nothing
        resultfilepath = joinpath(outpath, resultfilename)
        try
            if isfile(resultfilepath)
                open(resultfilepath) do file
                    result = JSON.parse(file)
                end
            end
        catch y
            result = nothing
        end

        if result != nothing
            ρ = [ float(x["re"]) + 1im * float(x["im"]) for x in result["rho"] ]
            t = [ float(x["re"]) + 1im * float(x["im"]) for x in result["t"] ]
            Γ = [ float(x["re"]) + 1im * float(x["im"]) for x in result["Gamma"] ]
            Δ = [ float(x["re"]) + 1im * float(x["im"]) for x in result["Delta"] ]
            return (result["BatchRun"], HFBSolution(ρ,t,Γ,Δ))
        end
    end

    return (0, nothing)
end

function runloop(solver ::HFBSolver;
                 outpath ::AbstractString="out",
                 nwarmup ::Integer = 200,
                 nbunch  ::Integer = 100,
                 nbatch  ::Integer = 500,
                 tolerance::Float64 = sqrt(eps(Float64)),
                 noise::Real = 0,
                 verbose::Bool=false,
                 update::Function=simpleupdate,
                 getnextsolution::Function=getnextsolution,
                 loop::Function=loop)

    @info "Entering runloop"
    progressbar = verbose && haskey(ENV, "TERM") && (ENV["TERM"] != "dumb")

    @assert(nwarmup >= 0)
    @assert(nbunch > 0)
    @assert(nbatch >= 0)

    @info "Making path $(outpath)"
    mkpath(outpath)

    hamspec = solver.hamiltonian
    uc = hamspec.unitcell

    @info "Making solution"
    currentsolution = newhfbsolution(solver.hfbcomputer)

    @info "Detecting result"
    batchrun_offset, lastsolution = detectresult(outpath)

    if lastsolution != nothing
        if iscompatible(currentsolution, lastsolution)
            @info("Using previous solution.")
            currentsolution = lastsolution
            @info ("Successfully read.")

            if verbose
                !isempty(currentsolution.ρ) && @info("  max |rho|  : $(maximum(abs.(currentsolution.ρ)))")
                !isempty(currentsolution.t) && @info("  max |t|    : $(maximum(abs.(currentsolution.t)))")
                !isempty(currentsolution.Γ) && @info("  max |Gamma|: $(maximum(abs.(currentsolution.Γ)))")
                !isempty(currentsolution.Δ) && @info("  max |Delta|: $(maximum(abs.(currentsolution.Δ)))")
            end

            if abs(noise) > 0
                @info("Adding noise to the values")
                currentsolution.ρ[:] += noise * (2*rand(Float64, size(currentsolution.ρ)) - 1)
                currentsolution.t[:] += noise * (2*rand(Complex128, size(currentsolution.t)) - 1 - 1im)
                currentsolution.Γ[:] += noise * (2*rand(Float64, size(currentsolution.Γ)) - 1)
                currentsolution.Δ[:] += noise * (2*rand(Complex128, size(currentsolution.Δ)) - 1 - 1im)
            end
        else
            @warn("Previous solution has wrong size. Using new random solution.")
            batchrun_offset = 0
            randomize!(solver.hfbcomputer, currentsolution)
        end
    else
        @info "Result not found"
        batchrun_offset = 0
        randomize!(solver.hfbcomputer, currentsolution)
    end

    if batchrun_offset == 0
        @info("Warming up")
        if progressbar
            @showprogress for run in 1:nwarmup
                currentsolution = getnextsolution(solver, currentsolution)
                currentsolution.Γ[:] = 0
            end
        else
            for run in 1:nwarmup
                currentsolution = getnextsolution(solver, currentsolution)
                currentsolution.Γ[:] = 0
            end
        end
    end

    for batchrun in 1:nbatch
        verbose && @info("BATCH $(batchrun + batchrun_offset)")
        callback = if progressbar
            p = Progress(nbunch)
            (i, n) -> next!(p)
        else
            (i, n) -> nothing
        end

        starttime = now()
        previoussolution = copy(currentsolution)
        currentsolution = loop(solver, currentsolution, nbunch;
                               update=update, callback=callback)
        endtime = now()
        verbose && @info("Duration: $(endtime - starttime)")

        maxdiff1 = isempty(currentsolution.ρ) ? 0.0 : maximum(abs.(currentsolution.ρ - previoussolution.ρ))
        maxdiff2 = isempty(currentsolution.t) ? 0.0 : maximum(abs.(currentsolution.t - previoussolution.t))
        maxdiff = max(maxdiff1, maxdiff2)

        (E, S, Ω) = hfbfreeenergy(solver, currentsolution; update=update)
        if verbose
            @info("Grand Potential = $Ω")
            @info("MaxDiff_rho = $maxdiff1")
            @info("MaxDiff_t   = $maxdiff2")
        end

        if isfile(joinpath(outpath, "result.json"))
            mv(joinpath(outpath, "result.json"),
               joinpath(outpath, "previous_result.json"),
               remove_destination=true)
        end
        open(joinpath(outpath, "result.json"), "w") do file
            outdict = OrderedDict([
                                   ("BatchRun",    batchrun_offset + batchrun),
                                   ("rho",         currentsolution.ρ),
                                   ("t",           currentsolution.t),
                                   ("Gamma",       currentsolution.Γ),
                                   ("Delta",       currentsolution.Δ),
                                   ("GrandPotential", Ω),
                                   ("MaxDiff_rho", maxdiff1),
                                   ("MaxDiff_t",   maxdiff2),
                                   ("StartTime",   Dates.format(starttime, "yyyy-mm-ddTHH:MM:SS.s")),
                                   ("EndTime",     Dates.format(endtime, "yyyy-mm-ddTHH:MM:SS.s")),
                                  ])
            JSON.print(file, outdict)
        end

        if maxdiff < tolerance
            return
        end

    end
end
