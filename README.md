# HartreeFockBogoliubov.jl

[![Build Status][travis-img]][travis-url]
[![Code Coverage][codecov-img]][codecov-url]

## Overview

HartreeFockBogoliubov.jl is a [Julia](https://julialang.org) package for solving interacting fermion Hamiltonian using Hartree-Fock-Bogoliubov (HFB) approach. This project aims to assist fast development of HFB solver for a generic Hamiltonian of the form
```math
\mathcal{H} = \sum_{\alpha,\beta} c_{\alpha}^{\dagger} t_{\alpha\beta} c_{\beta}
+ \frac{1}{2} \sum_{\alpha,\beta,\gamma,\delta} V_{\alpha\beta\gamma\delta} c_{\alpha}^{\dagger} c_{\beta}^{\dagger} c_{\delta} c_{\gamma}
```

## Installation

HartreeFockBogoliubov.jl is currently not registered in the [Julia Package Listing](https://pkg.julialang.org). To install, type
```
add https://github.com/kyungminlee/HartreeFockBogoliubov.jl
```
to Julia's package manager.


## Documentation

[![**STABLE**][docs-stable-img]][docs-stable-url] [![**LATEST**][docs-latest-img]][docs-latest-url]

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://kyungminlee.org/HartreeFockBogoliubov.jl/stable
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://kyungminlee.org/HartreeFockBogoliubov.jl/latest

[travis-img]: https://travis-ci.org/kyungminlee/HartreeFockBogoliubov.jl.svg?branch=master
[travis-url]: https://travis-ci.org/kyungminlee/HartreeFockBogoliubov.jl

[codecov-img]: https://codecov.io/gh/kyungminlee/HartreeFockBogoliubov.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/kyungminlee/HartreeFockBogoliubov.jl
