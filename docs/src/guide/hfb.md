# HFB: Hartree-Fock-Bogoliubov

## HFBComputer

`HFBComputer` is a struct which is used to calculate the Hartree-Fock-Bogoliubov Hamiltonian from mean-field parameters, and the parameters from the eigenstates of Hartree-Fock-Bogoliubov Hamiltonian.

```julia
mutable struct HFBComputer{O}
    unitcell ::UnitCell{O}
    hoppings ::Vector{Spec.Hopping}
    temperature ::Float64

    fermi ::Function
    ρ_registry ::Vector{CollectRow}
    t_registry ::Vector{CollectRow}
    Γ_registry ::Vector{DeployRow}
    Δ_registry ::Vector{DeployRow}
end
```

The type `CollectRow` keeps record of the "observables" that need to be measured in order to calculate the mean-field parameters. They are the ρ(i,j) and t(i,j) as defined by Goodman.
```julia
const CollectRow = Tuple{Bool, Int, Int, Vector{Float64}}
```

`DeployRow` contains necessary information to calculate Γ(i,j) and Δ(i,j), using the measured ρ(i,j) and t(i,j)
```julia
const DeployRow = Tuple{Bool, Int, Int, Vector{Float64}, Vector{Tuple{Int, Complex{Float64}, Bool}}}
```


## HFBSolver

`HFBSolver` is blah.

```julia
mutable struct HFBSolver{O}
    # Originals
    hamiltonian ::FullHamiltonian{O}
    size ::Vector{Int}
    temperature ::Float64

    # Derivatives
    momentumgrid ::Array{Vector{Float64}}
    hfbhamiltonian ::HFBHamiltonian{O}
    hfbcomputer ::HFBComputer{O}
    greencollectors ::Function
end
```

HFBSolver contains:

1. The full interacting Hamiltonian,
2. Size of the system (i.e. number of k-points in the Brillouin zone)
3. Temperature (needed to calculate the ρ and t)

And using these properties,
BLAH.


### Using HFBSolver to Find Self-Consistent Solution

```julia
loop(solver, solution, 100)
```
