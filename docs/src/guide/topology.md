# Topology

## Chern Number

To compute the total Chern number of a set of bands
```julia
unitcell = make_unitcell(...)
hoppings = [...]
n1 = 256
n2 = 256
selectbands = [1,2,3]
ch = chernnumber(unitcell, hoppings, n1, n2, selectband; tol=tol)
```
which breaks up the two-dimensional

## Z2 Index

Z2 index of time-reversal invariant topological insulator/superconductor.

```julia
z2 = z2index(...)
```
