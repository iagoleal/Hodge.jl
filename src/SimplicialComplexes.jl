# Topological / combinatorial operations on simplicial complexes.
# The functions in here do not see the recursive node structure
# but operate on a higher level of abstraction than those on the SimplexTree module.
# In special, in here we are always working with an entire simplicial complex
# and not with single nodes of a simplex tree
import .SimplexTrees
import LinearAlgebra

using Base: insert!

"""
    SimplicialComplex

Represent an abstract simplicial complex
where all vertices are labeled by integers.

Internally, use the [`SimplexTrees.SimplexTree`](@ref)
data structure.
"""
mutable struct SimplicialComplex
    dimension :: Int
    simplices :: SimplexTrees.SimplexTree

    SimplicialComplex() = new(-1, SimplexTrees.STRoot())
end


################################
# Construct simplicial complex #
################################

# Construct a simplicial complex from an iterator of simplices
function SimplicialComplex(iter)
    sc = SimplicialComplex()
    foreach(s -> insert!(sc, s), iter)
    return sc
end

Base.copy(sc::SimplicialComplex) = SimplicialComplex(simplices(sc))

"""
    insert!(sc::SimplicialComplex, σ)

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

"""
    skeleton(sc, k)

Return the k-skeleton of a [`SimplicialComplex`](@ref).
That is,
a SimplicialComplex with the same simplices as `sc`
up to dimension `k`.
"""
function skeleton(sc::SimplicialComplex, k::Integer)
    return SimplicialComplex(Iterators.flatten(simplices(sc,i) for i in 0:k))
end

#########################################
# Access Simplicial Complex information #
#########################################

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
function simplices end

function simplices(sc::SimplicialComplex) :: Vector{Vector{Int}}

        proper_faces = Iterators.flatten(SimplexTrees.getsimplices(sc.simplices, k) for k in 0:dimension(sc))
        return pushfirst!(collect(proper_faces), []) # Remember to add the empty face
end

function simplices(sc::SimplicialComplex, dim::Integer) :: Vector{Vector{Int}}
    if dim < 0 || dim > dimension(sc)
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
function numsimplices end

function numsimplices(sc::SimplicialComplex) :: Int
    return 1 + sum(Iterators.flatten(SimplexTrees.numsimplices(sc.simplices, k) for k in 0:dimension(sc)))
end

function numsimplices(sc::SimplicialComplex, dim::Integer) :: Int
    if dim < 0 || dim > dimension(sc)
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


########################
# Topological Features #
########################

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

"""
    betti(K)

Return the Betti numbers of a [`SimplicialComplex`](@ref).
"""
function betti end;

function betti(sc::SimplicialComplex, k::Integer) :: Int
    if k < 0 || k > dimension(sc)
        return 0
    end
    A = matrixify(laplacian, sc, k)
    return numsimplices(sc, k) - LinearAlgebra.rank(A)
end

betti(sc::SimplicialComplex) = [betti(sc,k) for k in 0:dimension(sc)]
