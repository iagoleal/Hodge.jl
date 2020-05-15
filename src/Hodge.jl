module Hodge

# Data structure for manipulation of simplicial complexes
include("SimplexTrees.jl")

# Data strucute for manipulation of sparse tensors
include("SparseTensors.jl")

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
