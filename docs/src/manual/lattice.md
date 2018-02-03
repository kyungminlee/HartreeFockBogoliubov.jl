# Lattice

## CarteCoord and FractCoord

## UnitCell

```julia
unitcell = newunitcell([1.0 0.0; 0.0 1.0]; OrbitalType=Tuple{Symbol, Symbol})
addorbital!(unitcell, (:A, :UP), FractCoord([0,0], [0.0, 0.0]))
addorbital!(unitcell, (:A, :DN), FractCoord([0,0], [0.0, 0.0]))
```
