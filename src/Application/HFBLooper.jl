using JSON
using DataStructures
using ArgParse
import YAML

using ProgressMeter
using HartreeFockBogoliubov
import HartreeFockBogoliubov.Spec
import HartreeFockBogoliubov.Generator
using HartreeFockBogoliubov.HFB


function detectresult(outpath)

    for resultfilename in ["result_concise.json", "result_history.json"]
        result = nothing
        resultfilepath = joinpath(outpath, resultfilename)
        if isfile(resultfilepath)
            open(resultfilepath) do file
                result = JSON.load(file)
            end
        end
        if result != nothing
            ρ = [ float(x["re"]) + 1im * float(x["im"]) for x in result["rho"] ]
            t = [ float(x["re"]) + 1im * float(x["im"]) for x in result["t"] ]
            Γ = [ float(x["re"]) + 1im * float(x["im"]) for x in result["Gamma"] ]
            Δ = [ float(x["re"]) + 1im * float(x["im"]) for x in result["Delta"] ]
            return HFBSolution(ρ,t,Γ,Δ)
        end
    end

    for resultfilename in ["result_concise.yaml", "result_history.yaml"]
        result = nothing
        resultfilepath = joinpath(outpath, resultfilename)
        if isfile(resultfilepath)
            open(resultfilepath) do file
                prev_yaml = YAML.load_all(file)
                for s in prev_yaml
                    result = s
                end
            end
        end

        if result != nothing
            ρ = [ float(x["real"]) + 1im * float(x["imag"]) for x in result["rho"] ]
            t = [ float(x["real"]) + 1im * float(x["imag"]) for x in result["t"] ]
            Γ = [ float(x["real"]) + 1im * float(x["imag"]) for x in result["Gamma"] ]
            Δ = [ float(x["real"]) + 1im * float(x["imag"]) for x in result["Delta"] ]
            return HFBSolution(ρ,t,Γ,Δ)
        end
    end
    return nothing
end

function runloop(solver ::HFBSolver;
    outpath ::AbstractString="out",
    nwarmup ::Integer = 200,
    nbunch  ::Integer = 100,
    nbatch  ::Integer = 500)
    assert(nwarmup >= 0)
    assert(nbunch > 0)
    assert(nbatch >= 0)
    mkpath(outpath)

    hamspec = solver.hamiltonian
    uc = hamspec.unitcell

    currentsolution = newhfbsolution(solver.hfbcomputer)
    lastsolution = detectresult(outpath)
    if lastsolution != nothing
        if iscompatible(currentsolution, lastsolution)
            println("Previous Solution matches in size with current solution")
            println("Let's use it")

            currentsolution = lastsolution

            println("Successfully read.")
            println("  maximum |rho|  : ", maximum(abs.(currentsolution.ρ)))
            println("  maximum |t|    : ", maximum(abs.(currentsolution.t)))
            println("  maximum |Gamma|: ", maximum(abs.(currentsolution.Γ)))
            println("  maximum |Delta|: ", maximum(abs.(currentsolution.Δ)))
        else
            println("Previous solution has different size from current solution")
            println("Something is wrong. Just using the new solution")
        end
    else
        randomize!(solver.hfbcomputer, currentsolution)
    end

    function nogammaupdate(sol::HFBSolution, newsol::HFBSolution)
        simpleupdate(sol, newsol)
        sol.Γ[:] = 0
        sol
    end

    println("Warm Up")
    for run in 1:nwarmup
        currentsolution = getnextsolution(solver, currentsolution)
        currentsolution.Γ[:] = 0
    end


    for batchrun in 1:nbatch
        println("BATCH $batchrun")
        nrun = nbunch + (batchrun % 2);
        p = Progress(nrun)

        callback = verbose ? (i, n) -> next!(p) : (i, n) -> nothing

        starttime = now()
        previoussolution = copy(currentsolution)
        currentsolution = loop(solver,
                               currentsolution,
                               nbunch + (batchrun % 2);
                               update=nogammaupdate,
                               callback=callback,
                               )
        endtime = now()
        println("Duration: ", (endtime - starttime))

        maxdiff1 = maximum(abs.(currentsolution.ρ - previoussolution.ρ))
        maxdiff2 = maximum(abs.(currentsolution.t - previoussolution.t))
        maxdiff = max(maxdiff1, maxdiff2)

        println("MaxDiff_rho = $maxdiff1")
        println("MaxDiff_t   = $maxdiff2")
        println("MaxDiff     = $maxdiff")

        open(joinpath(outpath, "result_history.yaml"), "a") do file
            FMT(x...) = begin
                foreach(z -> mydump(file, z), x)
                println(file)
            end
            FMT("---")
            FMT("BatchRun: ",    batchrun)
            FMT("rho: ",         currentsolution.ρ)
            FMT("t: ",           currentsolution.t)
            FMT("Gamma: ",       currentsolution.Γ)
            FMT("Delta: ",       currentsolution.Δ)
            FMT("MaxDiff_rho: ", maxdiff1)
            FMT("MaxDiff_t: ",   maxdiff2)
            FMT("StartTime: ",   starttime)
            FMT("EndTime: ",     endtime)
            FMT("...")
        end

        open(joinpath(outpath, "result_concise.yaml"), "w") do file
            FMT(x...) = begin
                foreach(z -> mydump(file, z), x)
                println(file)
            end
            FMT("---")
            FMT("BatchRun: ",    batchrun)
            FMT("rho: ",         currentsolution.ρ)
            FMT("t: ",           currentsolution.t)
            FMT("Gamma: ",       currentsolution.Γ)
            FMT("Delta: ",       currentsolution.Δ)
            FMT("MaxDiff_rho: ", maxdiff1)
            FMT("MaxDiff_t: ",   maxdiff2)
            FMT("StartTime: ",   starttime)
            FMT("EndTime: ",     endtime)
            FMT("...")
        end
    end
end
