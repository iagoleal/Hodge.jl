# Reference

## Simplicial Complexes

### Construction

```@docs
SimplicialComplex
insert!(::SimplicialComplex, s)
skeleton
```

### Accessing

```@docs
dimension
hassimplex
vertices
numvertices
simplices
numsimplices
```

### Topological Operators

```@docs
euler_characteristic
betti
```

## Cochains

### Construction

```@docs
Cochain
zero_cochain
indicator_cochain
```

### Operators

```@docs
basespace
basering
degree
norm
norm2
inner
cup
coboundary
coboundary_adj
laplacian
hodge
```
