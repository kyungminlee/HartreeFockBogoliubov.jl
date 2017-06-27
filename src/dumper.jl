export dumpall
export mydump


function mydump(io::IO, x::Integer)
  @printf(io, "%d", x)
end

function mydump(io::IO, x::AbstractFloat)
  @printf(io, "%a", x)
end

function mydump(io::IO, x::Complex{T}) where T <:Real
  print(io, "{ real: ")
  mydump(io, real(x))
  print(io, ", imag: ")
  mydump(io, imag(x))
  print(io, "}")
end

function mydump(io::IO, arr::AbstractVector{T}) where T
  print(io, "[")
  n = length(arr)
  first = true
  for i = 1:n
    if !first
      print(io, ", ")
    end
    mydump(io, arr[i])
    first = false
  end
  print(io, "]")
end

function mydump(io::IO, tup::Tuple)
  arr = collect(tup)
  mydump(io, arr)
end

function mydump(io::IO, str::AbstractString)
  print(io, str)
end

function mydump(io ::IO, x::Symbol)
  print(io, string(x))
end

function mydump(io::IO, x::UnitCell{T}) where T
  FMT(x...) = begin
    foreach(z -> mydump(io, z), x)
    println(io)
  end
  dim = dimension(uc)
  #=
  Dict( "Dimension" => dim,
        "LatticeVectors" => [uc.latticevectors[:, i] for i in 1:dim],
        "OrbitalType" => string(T),
        "Orbitals" => [
          Dict("Index" => idx, "Name" => name, "Coord" => )
        ]

        )
  =#

  FMT(    "\"UnitCell\": {")
  FMT(    "  \"Dimension\": $(dim),")
  FMT(    "  \"LatticeVectors\": [")
  for i=1:dim
    FMT(  "    ", uc.latticevectors[:, i])
  end
  FMT(    "  OrbitalType : \"$(T)\"")
  FMT(    "  Orbitals:")
  for (idx, (name, coord)) in enumerate(uc.orbitals)
    FMT(  "    - Index: $idx")
    FMT(  "      Name: ", name)
    FMT(  "      Coord: { whole: ", coord.whole,
                          ", fraction: ", coord.fraction, " }")
  end
end


function dumpall(io::IO, uc::UnitCell{T}) where T
  FMT(x...) = begin
    foreach(z -> mydump(io, z), x)
    println(io)
  end
  dim = dimension(uc)

  FMT(    "UnitCell:")
  FMT(    "  Dimension: $(dim)")
  FMT(    "  LatticeVectors:")
  for i=1:dim
    FMT(  "    - ", uc.latticevectors[:, i])
  end
  FMT(    "  OrbitalType : \"$(T)\"")
  FMT(    "  Orbitals:")
  for (idx, (name, coord)) in enumerate(uc.orbitals)
    FMT(  "    - Index: $idx")
    FMT(  "      Name: ", name)
    FMT(  "      Coord: { whole: ", coord.whole,
                          ", fraction: ", coord.fraction, " }")
  end
end

function dumpall(io::IO, solver::HFB.HFBSolver{T}) where T
  FMT(x...) = begin
    foreach(z -> mydump(io, z), x)
    println(io)
  end
  FMT(    "Registry:")
  FMT(    "  rho:")
  for (idx, (i, j, r)) in enumerate(solver.hfbcomputer.ρ_registry)
    FMT(  "    - Index: $idx")
    FMT(  "      Row: $i")
    FMT(  "      Col: $j")
    FMT(  "      Vec: ", r)
  end
  FMT(    "  t:")
  for (idx, (i, j, r)) in enumerate(solver.hfbcomputer.t_registry)
    FMT(  "    - Index: $idx")
    FMT(  "      Row: $i")
    FMT(  "      Col: $j")
    FMT(  "      Vec: ", r)
  end
  FMT(    "  Gamma:")
  for (idx, (i, j, r, srcs)) in enumerate(solver.hfbcomputer.Γ_registry)
    FMT(  "    - Index: $idx")
    FMT(  "      Row: $i")
    FMT(  "      Col: $j")
    FMT(  "      Vec: ", r)
    FMT(  "      Sources:")
    for (srcidx, amplitude, cj) in srcs
      FMT("        - Index: $srcidx")
      FMT("          Amplitude: ", amplitude)
      FMT("          Conjugate: $cj")
    end
  end
  FMT(    "  Delta:")
  for (idx, (i, j, r, srcs)) in enumerate(solver.hfbcomputer.Δ_registry)
    FMT(  "    - Index: $idx")
    FMT(  "      Row: $i")
    FMT(  "      Col: $j")
    FMT(  "      Vec: ", r)
    FMT(  "      Sources:")
    for (srcidx, amplitude, minussign) in srcs
      FMT("        - Index: $srcidx")
      FMT("          Amplitude: ", amplitude)
      FMT("          Sign: $minussign")
    end
  end
  kpoints = reshape(solver.momentumgrid, length(solver.momentumgrid))
  FMT(    "  Momenta:")
  for k in kpoints
    FMT(  "    - ", k)
  end
end

function dumpall(io::IO,
                 hamspec ::Spec.SpecHamiltonian{T},
                 solver ::HFB.HFBSolver{T},
                 currentsolution ::HFB.HFBSolution,
                 previoussolution ::HFB.HFBSolution) where T
  uc = hamspec.unitcell

  FMT(x...) = begin
    foreach(z -> mydump(io, z), x)
    println(io)
  end

  dim = dimension(uc)
  dumpall(io, uc)
  dumpall(io, solver)

  FMT(    "CurrentSolution:")
  FMT(    "  rho:")
  for ρ in currentsolution.ρ
    FMT(  "    - ", ρ)
  end
  FMT(    "  t:")
  for t in currentsolution.t
    FMT(  "    - ", t)
  end
  FMT(    "  Gamma:")
  for Γ in currentsolution.Γ
    FMT(  "    - ", Γ)
  end
  FMT(    "  Delta:")
  for Δ in currentsolution.Δ
    FMT(  "    - ", Δ)
  end

  FMT(    "PreviousSolution:")
  FMT(    "  rho:")
  for ρ in previoussolution.ρ
    FMT(  "    - ", ρ)
  end
  FMT(    "  t:")
  for t in previoussolution.t
    FMT(  "    - ", t)
  end
  FMT(    "  Gamma:")
  for Γ in previoussolution.Γ
    FMT(  "    - ", Γ)
  end
  FMT(    "  Delta:")
  for Δ in previoussolution.Δ
    FMT(  "    - ", Δ)
  end
end
