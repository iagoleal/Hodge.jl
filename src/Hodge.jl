module Hodge

# Data structures for manipulation of simplicial complexes
include("SimplexTrees.jl")

# Topological / Geometrical structures
include("SimplicialComplexes.jl")

export SimplicialComplex
export dimension,
       hassimplex,
       vertices,
       numvertices,
       simplices,
       numsimplices,
       insert!

end # module
