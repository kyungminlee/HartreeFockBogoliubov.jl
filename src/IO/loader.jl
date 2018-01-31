import YAML

function myload(text)
    data = YAML.load(text)
    data = myconv(data)
    return data
end

function myconv(obj)
    const hexpat = r"[+-]?0x([0-9a-fA-F]+(\.[0-9a-f]*)?|(\.[0-9a-f]*))(p[+-]?[0-9a-f]+)?"i
    const complexkeyword = sort(["real", "imag"])
    if isa(obj, Dict)
        #if sort(keys(obj)) == complexkeyword
        if haskey(obj, "real") && haskey(obj, "imag")
            return complex(myconv(obj["real"]), myconv(obj["imag"]))
        else
            return Dict(k => myconv(v) for (k,v) in obj)
        end
    elseif isa(obj, AbstractArray)
        return [myconv(v) for v in obj]
    elseif isa(obj, Integer)
        return obj
    elseif isa(obj, AbstractFloat)
        return obj
    elseif isa(obj, AbstractString)
        if ismatch(hexpat, obj)
            return parse(Float64, obj)
        else
            return obj
        end
    else
        return obj
    end
end

tuplify(obj::Vector) = (map(tuplify, obj)...)
tuplify(obj::Tuple) = map(tuplify, obj)
tuplify(obj) = obj


function load_parameter_file(data)
    data["Parameters"]
    unitcell = data["UnitCell"]
    ρ_registry = load_ρ_registry(data["Registry"]["rho"])
    t_registry = load_t_registry(data["Registry"]["t"])
    Γ_registry = load_ρ_registry(data["Registry"]["Gamma"])
    Δ_registry = load_ρ_registry(data["Registry"]["Delta"])


end


function load_unitcell(dict)
    dim ::Int= dict["Dimension"]
    if dim < 0
        error("Dimension should be non-negative")
    end
    lattice_vectors ::Matrix{Float64} = dict["LatticeVectors"]

    orbital_type = eval(parse(dict["OrbitalType"]))

    isa(orbital_type, Type) || error("orbital type should be a valid Julia type")

    unitcell = newunitcell(lattice_vectors; OrbitalType=orbital_type)

    orbitals ::Vector{Tuple{orbital_type, FractCoord}} = []
    for (index2, orbital) in enumerate(dict["Orbitals"])
        index ::Int = orbital["Index"]
        name = orbital["Name"]
        coord = orbital["Coord"]

        @assert(index == index2)

        if orbital_type <: Tuple
            name = convert(orbital_type, tuplify(name))
        end
        fractcoord = FractCoord(coord["whole"], coord["fraction"])
        addorbital!(unitcell, name, fractcoord)
    end
    return unitcell
end

function load_ρ_registry(data)
    ρ_registry = CollectRow[]
    for (index2, item) in enumerate(data)
        index ::Int = item["Index"]
        diag ::Bool = item["Diag"]
        row ::Int = item["Row"]
        col ::Int = item["Col"]
        vec ::Vector{Float64} = item["Vec"]
        @assert(index == index2)
        push!(ρ_registry, CollectRow(diag, row, col, vec))
    end
    return ρ_registry
end

function load_t_registry(data)
    t_registry = CollectRow[]
    for (index2, item) in enumerate(data)
        index ::Int = item["Index"]
        diag ::Bool = item["Diag"]
        row ::Int = item["Row"]
        col ::Int = item["Col"]
        vec ::Vector{Float64} = item["Vec"]
        @assert(index == index2)
        @assert(diag == false)
        push!(t_registry, CollectRow(diag, row, col, vec))
    end
    return t_registry
end

function load_Γ_registry(data)
    Γ_registry = DeployRow[]
    for (index2, item) in enumerate(data)
        index ::Int = item["Index"]
        diag ::Bool = item["Diag"]
        row ::Int = item["Row"]
        col ::Int = item["Col"]
        vec ::Vector{Float64} = item["Vec"]
        @assert(index == index2)
        sources = Tuple{Int64, Complex128, Bool}[
        (src["Index"], src["Amplitude"], src["Conjugate"])
        for src in item["Sources"]
            ]
            push!(Γ_registry, DeployRow(diag, row, col, vec, sources))
        end
        return Γ_registry
    end

    function load_Δ_registry(data)
        Δ_registry = DeployRow[]
        for (index2, item) in enumerate(data)
            index ::Int = item["Index"]
            diag ::Bool = item["Diag"]
            row ::Int = item["Row"]
            col ::Int = item["Col"]
            vec ::Vector{Float64} = item["Vec"]
            @assert(index == index2)
            sources = Tuple{Int64, Complex128, Bool}[
            (src["Index"], src["Amplitude"], src["Conjugate"])
            for src in item["Sources"]
                ]
                push!(Δ_registry, DeployRow(diag, row, col, vec, sources))
            end
            return Δ_registry
        end
