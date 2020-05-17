# Topological / combinatorial operations on simplicial complexes.
# The functions in here do not see the recursive node structure
# but operate on a higher level of abstraction than those on the SimplexTree module.
# In special, in here we are always working with an entire simplicial complex
# and not with single nodes of a simplex tree
import .SimplexTrees

using Base: insert!

"""
    SimplicialComplex

Represent an abstract simplicial complex
where all vertices are labeled by integers.

Internally, use the [`SimplexTree`](@ref)
data structure.
"""
mutable struct SimplicialComplex
    dimension :: Int
    simplices :: SimplexTrees.SimplexTree

    SimplicialComplex() = new(-1, SimplexTrees.STRoot())
end

# Construct a simplicial complex from an iterator of simplices
function SimplicialComplex(iter)
    sc = SimplicialComplex()
    foreach(s -> insert!(sc, s), iter)
    return sc
end

# Acessor functions
"""
    dimension(sc)

Return the dimension of the largest simplex
on a [`SimplicialComplex`](@ref).
"""
dimension(sc::SimplicialComplex) = sc.dimension

# Access simplices

"""
    vertices(sc)

Return an array containg the vertices of the [`SimplicialComplex`](@ref) `sc`.

See also [`simplices`](@ref).
Notice that this function returns an array of integers
while `simplices(sc, 0)` returns an array of singleton arrays
containing the vertices.
"""
@inline function vertices(sc::SimplicialComplex)
    return SimplexTrees.children_label(sc.simplices)
end

"""
    numvertices(sc)

Return how many vertices (1-simplices)
the [`SimplicialComplex`](@ref) `sc` has.
"""
@inline function numvertices(sc::SimplicialComplex)
    return length(sc.simplices.children)
end

"""
    simplices(sc[, k])

Return all simplices of the [`SimplicialComplex`](@ref)
`sc` whose dimension equals `k`.

If the parameter `k` is not given,
return __all__ simplices of `sc`
including the empty face.
"""
function simplices(sc::SimplicialComplex, dim=nothing)
    if dim === nothing
        proper_faces = vcat(map( k -> SimplexTrees.getsimplices(sc.simplices, k)
                       , 0:dimension(sc))...)
        return pushfirst!(proper_faces, []) # Remember to add the empty face
    elseif dim < 0 || dim > dimension(sc)
        return Vector{Int}[]
    else
        return SimplexTrees.getsimplices(sc.simplices, dim)
    end
end

"""
    numsimplices(sc[, k])

Return the number of k-dimensional simplices of `sc`.

If the parameter `k` is not given,
return the total number of simplices
including the empty face.

This function is a more efficient implementation of
`length ∘ simplices`.
"""
function numsimplices(sc::SimplicialComplex, dim=nothing)
    if dim === nothing
        return 1 + sum(map( k -> SimplexTrees.numsimplices(sc.simplices, k)
                          , 0:dimension(sc)))
    elseif dim < 0 || dim > dimension(sc)
        return 0
    elseif dim == 0
        return numvertices(sc)
    else
        return SimplexTrees.numsimplices(sc.simplices, dim)
    end
end

"""
    hassimplex(sc, σ)

Test whether the [`SimplicialComplex`](@ref) `sc`
contains the simplex `σ`.
"""
@inline function hassimplex(sc::SimplicialComplex, simplex)
    return SimplexTrees.hassimplex(sc.simplices, simplex)
end

## Construct simplicial complex

"""
    insert!(sc, σ)

Insert the simplex `σ` and all its faces
on the [`SimplicialComplex`](@ref) `sc`.

The simplex does not need to be ordered.
"""
function Base.insert!(sc::SimplicialComplex, simplex)
    SimplexTrees.insert!(sc.simplices, simplex)
    newdim = length(simplex) -1
    if newdim > sc.dimension
        sc.dimension = newdim
    end
    return sc
end

## Topological features
"""
    euler_characteristic(K)

Return the Euler characteristic of a [`SimplicialComplex`](@ref).

"""
function euler_characteristic(sc::SimplicialComplex)
    x  = 0 :: Int
    sg = 1 :: Int
    for i in 0:dimension(sc)
        x += sg * numsimplices(sc, i)
        sg = -sg
    end
    return x
end
