# HartreeFockBogoliubov.jl Documentation

## Overview

HartreeFockBogoliubov.jl is a [Julia](https://julialang.org) package for solving interacting fermion Hamiltonian using Hartree-Fock-Bogoliubov (HFB) approach. This project aims to assist fast development of HFB solver for a generic Hamiltonian of the form
```math
\mathcal{H} = \sum_{\alpha,\beta} c_{\alpha}^{\dagger} t_{\alpha\beta} c_{\beta}
+ \frac{1}{2} \sum_{\alpha,\beta,\gamma,\delta} V_{\alpha\beta\gamma\delta} c_{\alpha}^{\dagger} c_{\beta}^{\dagger} c_{\delta} c_{\gamma}
```

## Installation

`HartreeFockBogoliubov.jl` is not a registered package. Install it by

HartreeFockBogoliubov.jl is currently not registered in the [Julia Package Listing](https://pkg.julialang.org). To install, type
```julia
add https://github.com/kyungminlee/HartreeFockBogoliubov.jl
```
to Julia's package manager.

## Usage

To construct HFB solver
1. Specify a unitcell, i.e. the lattice constants.
2. Specify the "orbitals" of the unitcell, and their locations. Here the word "orbital" refers to any fermionic degrees of freedom, and includes sites, orbitals in the conventional sense, and spins, etc.
3. Specify the hoppings between the orbitals, together with their real-space displacements. The displacement is required since, with the periodic boundary condition, hopping from orbital α to orbital β within the same unitcell is different from the hopping from orbital α to β in different unitcells.
4. Specify the interactions. As for the hoppings, you need to specify the orbitals and the displacement of the interactions.
5. The full interacting Hamiltonian constructed so far can be passed to a `HFBSolver` which decomposes it into different HFB channels.

## License

```
MIT License

Copyright (c) 2016 Kyungmin Lee

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
