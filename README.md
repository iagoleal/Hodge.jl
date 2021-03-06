# Hodge.jl

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![Stable Docs][docs-stable-img]][docs-stable-url] [![Latest Docs][docs-dev-img]][docs-dev-url] | [![Build Status][build-img]][build-url]  [![Code Coverage][codecov-img]][codecov-url] |

Library for manipulation of simplicial complexes
and cochains defined over them.

The focus of this package is on the formalism of
**Discrete Exterior Calculus**
and its applications to statistical rankings
and extraction of topological features from simplicial triangulations.

## Installation
To install this package, all you have to do is to enter `]` on the Julia REPL and write
```julia
pkg> add Hodge
```

## Basic Usage
To use this package,
one starts defining a simplicial complex
```julia
using Hodge
K = SimplicialComplex([(1,2,3), (1,2,4), [1], [1,5,9,6], (2,6)])
```

Then it is possible to retrieve topological information from the complex
```julia
euler_characteristic(K)
betti(K)
```

Or one can define Cochains over `K`,
which are skew-symmetric tensors over the simplices of `K`,
and work with them
```julia
f = Cochain{Float64, 2}(K)
f[1,2,3] = 3.0
f[1,2,4] = -5.9
f[4,3,2] = 13.2
f[1,4,3] = π

g = coboundary(f)
w = coboundary_adj(f)
h = laplacian(f)
c = cup(f, h)

a, b, c = hodge(f)
```

See the [documentation][docs-stable-url] for a more comprehensive explanation.

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://iagoleal.github.io/Hodge.jl/stable/

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://iagoleal.github.io/Hodge.jl/dev/

[build-img]: https://github.com/iagoleal/Hodge.jl/actions/workflows/ci.yml/badge.svg?branch=master
[build-url]: https://github.com/iagoleal/Hodge.jl/actions/workflows/ci.yml

[codecov-img]: https://codecov.io/github/iagoleal/Hodge.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/iagoleal/Hodge.jl?branch=master
