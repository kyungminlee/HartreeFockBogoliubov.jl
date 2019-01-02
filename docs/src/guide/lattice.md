# Lattice

## CarteCoord and FractCoord

`CarteCoord` is simply an alias for `Vector{Float64}`, which represents the coordinates in the real-space.
`FractCoord` on the other hand, represents a location in the units of the lattice vectors.
For example, `fc` in the following code represents a location in two-dimensional space whose `CarteCoord` is
``2.1 \mathbf{a}_1 + 3.2 \mathbf{a}_2`` where ``\mathbf{a}_1`` and ``\mathbf{a}_2`` are lattice vectors.
```
fc = FractCoord([2, 3], [0.1, 0.2])
```
`FractCoord` can also be created in the following ways
```
fc2 = FractCoord([2.3, 3.2])
fc3 = FractCoord(2)
```
`fc2` is the same as `fc`, and `fc3` represents origin in two dimensions.


`FractCoord` supports the following arithmetic operations
```julia
fc = FractCoord([2, 3], [0.1, 0.2])

fc + fc
fc - fc
fc + [1, 0]
fc - [1, 0]
```


Conversion between `CarteCoord` and `FractCoord`:
```julia
latticevectors = [2.0 0.0; 0.0 1.0]

fc = FractCoord([2, 3], [0.1, 0.2])
cc = fract2carte(latticevectors, fc)
fc2 = carte2fract(latticevectors, cc)
```

`cc` is `[4.2, 3.2]` and `fc ≈ fc2` is `true`.


## UnitCell

A `UnitCell` is defined by lattice vectors and orbitals (orbitals here contain all possible degrees of freedom, such as spin, physical "orbital", sublattice, etc.).

```julia
unitcell = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=Tuple{Symbol, Symbol})
addorbital!(unitcell, (:A, :UP), FractCoord([0,0], [0.0, 0.0]))
addorbital!(unitcell, (:A, :DN), FractCoord([0,0], [0.0, 0.0]))
```
