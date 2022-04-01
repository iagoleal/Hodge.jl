# Hodge.jl

This package exports two main types,
[`SimplicialComplex`](@ref) and [`Cochain`](@ref),
together with methods to work with their topological and algebraic properties.

The topological operations on this package
are all done via the discrete Laplacian operator.
This includes the method [`betti`](@ref),
which calculates the Betti numbers of a simplicial complex,
and the method [`hodge`](@ref)
which calculates the discrete Hodge decomposition of a cochain.

## Installation
This package can be installed using the Julia Package Manager. Simply open the REPL, enter `]` and run

```julia
pkg> add Hodge
```

## Bibliography

`Hodge.jl` is based on a scientific initiation
that I did with [Prof. João Paixão](https://www.joaopaixao.com)
while an undergraduate at UFRJ.

The simplicial complex type is built upon the __Simplex Tree__ data structure,
described in the [paper](https://hal.inria.fr/hal-00707901v1/document):
- Jean-Daniel Boissonnat, Clément Maria. The Simplex Tree: An Efficient Data Structure for General Simplicial Complexes. [Research Report] RR-7993, 2012, pp.20. hal-00707901v1

The idea of representing cochains as skew-symmetric tensors
and using them to write the discrete Hodge decomposition was based on the
[paper](https://link.springer.com/article/10.1007%2Fs10107-010-0419-x):
- Jiang, X., Lim, L., Yao, Y. et al. Statistical ranking and combinatorial Hodge theory. Math. Program. 127, 203–244 (2011). https://doi.org/10.1007/s10107-010-0419-x
