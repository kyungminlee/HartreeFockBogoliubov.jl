module HFBLooper

export detectresult
export runloop

using Dates
using JSON
using DataStructures

using ProgressMeter
using HartreeFockBogoliubov
import HartreeFockBogoliubov.Spec
import HartreeFockBogoliubov.Generator
import HartreeFockBogoliubov.HFB
using HartreeFockBogoliubov.HFB


function detectresult(outpath ::AbstractString)
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
            return (result["BatchRun"], HFBAmplitude(ρ, t), HFBField(Γ, Δ))
        end
    end
    return (0, nothing, nothing)
end


function runloop(solver ::HFBSolver;
                 outpath ::AbstractString="out",
                 nwarmup ::Integer=200,
                 nbunch  ::Integer=100,
                 nbatch  ::Integer=500,
                 tolerance ::Real=sqrt(eps(Float64)),
                 noise ::Real=0,
                 verbose ::Bool=false,
                 update ::Function=simpleupdate,
                 loop ::Function=loop)

    verbose && @info "Entering runloop"
    progressbar = verbose && haskey(ENV, "TERM") && (ENV["TERM"] != "dumb")

    if nwarmup < 0
        throw(ArgumentError("nwarmup must be non-negative"))
    elseif nbunch <= 0
        throw(ArgumentError("nbunch must be positive"))
    elseif nbatch <= 0
        throw(ArgumentError("nbatch must be non-negative"))
    end

    verbose && @info "Making path $(outpath)"
    mkpath(outpath)

    hamspec = solver.hamiltonian
    uc = hamspec.unitcell

    verbose && @info "Making solution"
    current_hfbfield = make_hfbfield(solver)
    current_hfbamplitude = make_hfbamplitude(solver)

    verbose && @info "Detecting result"
    batchrun_offset, last_hfbamplitude, last_hfbfield = detectresult(outpath)

    if last_hfbfield != nothing
        if iscompatible(current_hfbamplitude, last_hfbamplitude) && iscompatible(current_hfbfield, last_hfbfield)
            verbose && @info "Using previous solution."
            current_hfbfield = last_hfbfield
            current_hfbamplitude = last_hfbamplitude
            verbose && @info "Successfully read."

            if verbose
                !isempty(current_hfbamplitude.ρ) && @info("  |ρ| : $(minimum(abs.(current_hfbamplitude.ρ))) -- $(maximum(abs.(current_hfbamplitude.ρ)))")
                !isempty(current_hfbamplitude.t) && @info("  |t| : $(minimum(abs.(current_hfbamplitude.t))) -- $(maximum(abs.(current_hfbamplitude.t)))")
                !isempty(current_hfbfield.Γ) && @info("  |Γ| : $(minimum(abs.(current_hfbfield.Γ))) -- $(maximum(abs.(current_hfbfield.Γ)))")
                !isempty(current_hfbfield.Δ) && @info("  |Δ| : $(minimum(abs.(current_hfbfield.Δ))) -- $(maximum(abs.(current_hfbfield.Δ)))")
            end

            if abs(noise) > 0
                verbose && @info("Adding noise to the values")
                current_hfbfield.Γ[:] += noise * (2*rand(Float64, size(current_hfbfield.Γ)) - 1)
                current_hfbfield.Δ[:] += noise * (2*rand(Complex128, size(current_hfbfield.Δ)) - 1 - 1im)
            end
        else
            @warn "Previous solution has wrong size. Using new random solution."
            batchrun_offset = 0
            randomize!(solver, current_hfbamplitude)
            current_hfbfield = make_hfbfield(solver, current_hfbamplitude)
        end
    else
        verbose && @info "Result not found"
        batchrun_offset = 0
        randomize!(solver, current_hfbamplitude)
        current_hfbfield = make_hfbfield(solver, current_hfbamplitude)
        fill!(current_hfbfield.Γ, 0)
    end

    callback = if progressbar
        p = Progress(nbunch)
        (i, n) -> next!(p)
    else
        (i, n) -> nothing
    end

    if batchrun_offset == 0 #TODO WHY?
        verbose && @info "Warming up"
        (current_hfbamplitude, current_hfbfield) = loop(solver, current_hfbfield, nwarmup;
                                                        update=update, callback=callback)
    end

    for batchrun in 1:nbatch
        verbose && @info("BATCH $(batchrun + batchrun_offset)")

        starttime = now()

        previous_hfbamplitude = deepcopy(current_hfbamplitude)
        previous_hfbfield = deepcopy(current_hfbfield)

        (current_hfbamplitude, current_hfbfield) = loop(solver, current_hfbfield, nbunch;
                                                        update=update, callback=callback)
        endtime = now()
        verbose && @info("Duration: $(endtime - starttime)")

        maxdiff_ρ = isempty(current_hfbamplitude.ρ) ? 0.0 : maximum(abs.(current_hfbamplitude.ρ - previous_hfbamplitude.ρ))
        maxdiff_t = isempty(current_hfbamplitude.t) ? 0.0 : maximum(abs.(current_hfbamplitude.t - previous_hfbamplitude.t))
        maxdiff_Γ = isempty(current_hfbfield.Γ) ? 0.0 : maximum(abs.(current_hfbfield.Γ - previous_hfbfield.Γ))
        maxdiff_Δ = isempty(current_hfbfield.Δ) ? 0.0 : maximum(abs.(current_hfbfield.Δ - previous_hfbfield.Δ))
        maxdiff = max(maxdiff_ρ, maxdiff_t)

        (E, S, Ω) = hfbfreeenergy(solver, current_hfbfield; update=update)
        if verbose
            max_ρ = isempty(current_hfbamplitude.ρ) ? 0.0 : maximum(abs.(current_hfbamplitude.ρ))
            max_t = isempty(current_hfbamplitude.t) ? 0.0 : maximum(abs.(current_hfbamplitude.t))
            max_Γ = isempty(current_hfbfield.Γ) ? 0.0 : maximum(abs.(current_hfbfield.Γ))
            max_Δ = isempty(current_hfbfield.Δ) ? 0.0 : maximum(abs.(current_hfbfield.Δ))
            @info("Grand potential = $Ω")
            @info("max(|Δρ|) = $maxdiff_ρ")
            @info("max(|Δt|) = $maxdiff_t")
            @info("max(|ΔΓ|) = $maxdiff_Γ")
            @info("max(|ΔΔ|) = $maxdiff_Δ")
            @info("max(|ρ|) = $max_ρ")
            @info("max(|t|) = $max_t")
            @info("max(|Γ|) = $max_Γ")
            @info("max(|Δ|) = $max_Δ")
        end

        if isfile(joinpath(outpath, "result.json"))
            mv(joinpath(outpath, "result.json"),
               joinpath(outpath, "previous_result.json"),
               force=true)
        end
        open(joinpath(outpath, "result.json"), "w") do file
            outdict = OrderedDict([
                                   ("BatchRun",    batchrun_offset+batchrun),
                                   ("rho",         current_hfbamplitude.ρ),
                                   ("t",           current_hfbamplitude.t),
                                   ("Gamma",       current_hfbfield.Γ),
                                   ("Delta",       current_hfbfield.Δ),
                                   ("GrandPotential", Ω),
                                   ("MaxDiff_rho", maxdiff_ρ),
                                   ("MaxDiff_t",   maxdiff_t),
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



end #module HFBLooper
