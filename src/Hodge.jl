module Hodge

# Data structure for manipulation of simplicial complexes
include("SimplexTrees.jl")

# Iterative method for solving linear equations
include("ConjugateGradient.jl")

# Topological / Geometrical structures
include("SimplicialComplexes.jl")

export SimplicialComplex
export dimension,
       hassimplex,
       vertices,
       numvertices,
       simplices,
       numsimplices,
       euler_characteristic
       # remove!,
       # flag,
       # betti,
       # expansion,
       # cofaces
       # skeleton

# Cochain algebra over a simplicial complex
include("Cochains.jl")

export Cochain
export basespace,
       degree,
       indicator_cochain,
       norm,
       norm2,
       inner,
       cup,
       coboundary,
       coboundary_adj,
       laplacian,
       hodge

end # module
