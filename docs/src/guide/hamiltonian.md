# Hamiltonian: Full Interacting Hamiltonian

## Terms of Hamiltonian

`HartreeFockBogoliubov.jl` currently supports Hamiltonian with hopping and (quartic) interaction terms.
The hopping terms can be grouped into "diagonal" terms and "offdiagonal" terms:
```math
\begin{align}
H_{\text{diagonal-hopping}} &= \sum_{i} t_{ii} c_{i}^{*} c_{i} \\
H_{\text{offdiagonal-hopping}} &= \sum_{i \neq j} t_{ij} c_{i}^{*} c_{j}
\end{align}
```
The hermiticity of the Hamiltonian requires that ``t_{ii}`` be real, and ``t_{ij} = t_{ji}^{*}``. Thus we can rewrite the offdiagonal term as
```math
H_{\text{offdiagonal-hopping}} = \sum_{i \lt j} t_{ij} c_{i}^{*} c_{j} + t_{ij}^{*} c_{j}^{*} c_{i}
```
To incorporate this we define two `struct`s: `HoppingDiagonal` and `HoppingOffdiagonal`:
``` julia
struct HoppingDiagonal{R<:Real}
    amplitude ::R
    i ::Int
    Ri ::Vector{Int}
end

struct HoppingOffdiagonal{C<:Number}
    amplitude ::C
    i ::Int
    j ::Int
    Ri ::Vector{Int}
    Rj ::Vector{Int}
end
```
`i` and `j` are integers representing the "index" of the orbital, and `Ri` and `Rj` represents which unitcell they are in.
The "unitcell coordinates" `Ri` and `Rj` are necessary for the representation of the Hamiltonian in the momentum space.
A single `HoppingOffdiagonal` represents both ``t_{ij} c_{i}^{*} c_{j}`` and its hermitian conjugate.


Similarly for the interaction, two structs are defined:
```julia
struct InteractionDiagonal{R<:Real}
    amplitude ::R
    i ::Int
    j ::Int
    Ri ::Vector{Int}
    Rj ::Vector{Int}
end

struct InteractionOffdiagonal{C<:Number}
    amplitude ::C
    i ::Int
    j ::Int
    k ::Int
    l ::Int
    Ri ::Vector{Int}
    Rj ::Vector{Int}
    Rk ::Vector{Int}
    Rl ::Vector{Int}
end
```
which represent
```math
U c_{i}^{*} c_{j}^{*} c_{j} c_{i}, \qquad i \lt j
```
and
```math
U     c_{i}^{*} c_{j}^{*} c_{l} c_{k} +
U^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}, \qquad
i \lt j \wedge k \lt l \wedge [ i \lt k \vee ( i = k \wedge j \lt l )]
```
In other words, (i,j) and (k,l) each need to be in lexicographical order, as well as
``(i,j) < (k,l)``.


The following unions are defined
```julia
const Hopping = Union{HoppingDiagonal, HoppingOffdiagonal}
const Interaction = Union{InteractionDiagonal, InteractionOffdiagonal}
```


### Functions Creating Terms

`HartreeFockBogoliubov.jl` provides functions that allows generation of hopping terms and interaction terms from the orbital names and their `CarteCoord`: `hoppingbycarte` and `interactionbycarte`.


### FullHamiltonian

The full interacting Hamiltonian is represented by the `FullHamiltonian` class
```julia
mutable struct FullHamiltonian{O}
    unitcell ::UnitCell{O}
    hoppings ::Vector{Hopping}
    interactions ::Vector{Interaction}
end
```

`addhopping!` and `addinteraction!`.
