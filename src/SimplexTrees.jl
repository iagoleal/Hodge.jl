## This module deals with the computational aspect
## of storing and manipulating a simplicial complex.
## All the geometric primitives are defined on SimplicialComplex.jl
## using the methods of this module.
module SimplexTrees

"""
Store a simplicial complex
using the Simplex Tree data structure.

Based on the paper:
 - Jean-Daniel Boissonnat, Clément Maria. The Simplex Tree: An Efficient Data Structure for General Simplicial Complexes. [Research Report] RR-7993, 2012, pp.20. hal-00707901v1
"""
abstract type SimplexTree end

# Node of the SimplexTree
# Represents a non-empty simplex
struct STNode <: SimplexTree
    label    :: Int
    children :: Array{STNode, 1}
    parent   :: SimplexTree
    STNode(label, parent) = new(label, [], parent)
end

# Tree root
# Represents empty face
struct STRoot <: SimplexTree
    children :: Array{STNode, 1}
    STRoot() = new([])
end

# Pertinence / propriety tests

"""
    isemptyface(node)

Test whether `node` represents an empty simplex.
"""
isemptyface(x::SimplexTree) = typeof(x) == STRoot

"""
    isleaf(node)

Test whether `node` is a leaf of a [`SimplexTree`](@ref)
"""
isleaf(x::SimplexTree) = isempty(x.children)

"""
    haschild(node, j)

Test whether `node` has a child with label `j`.
"""
function haschild(node::SimplexTree, val::Integer)
    return val in map(x -> x.label, node.children)
end

"""
    hassimplex(node, s)

Test whether the [`SimplexTree`](@ref) rooted at `node`
contains the simplex `s`.
"""
function hassimplex(complex::SimplexTree, simplex)
    if isempty(simplex)
        return true
    elseif !haschild(complex, simplex[1])
        return false
    else
        return hassimplex(child_from_label(complex, simplex[1]), simplex[2:end])
    end
end

## Getter functions
## Useful for converting between simplex and node information

"""
    children_label(node)

Return an array containing the labels of `node`'s children.
"""
function children_label(node::SimplexTree)
    return map(x -> x.label, node.children)
end

"""
    child_from_label(node, j)

Return the child of `node` whose label is `j`.
If no child has label `j`,
return `nothing`.
"""
function child_from_label(node::SimplexTree, val::Integer)
    label = findfirst(x -> x.label == val, node.children)
    if isnothing(label)
        return nothing
    else
        return node.children[label]
    end
end


"""
    getsimplices(st, k)

Return all the `k`-simplices of `complex`.
"""
function getsimplices(st::SimplexTree, dim::Integer)
    if dim < 0
        return Vector{Int}[]
    else
        return foldST_depth(st, dim+1;
                            basecase = x -> Vector{Int}[simplex_from_node(x)],
                            joiner   = x -> Vector{Vector{Int}}(vcat(x...))
                           )
    end
end

"""
    numsimplices(st, k)

Return the number of `k`-simplices of `st`.

This is a more efficient implementation of
`length ∘ getsimplices`
"""
function numsimplices(st::SimplexTree, dim::Integer)
    if dim < 0
        return 0
    else
        return foldST_depth(st, dim+1;
                            basecase = x -> 1,
                            joiner   = sum)
    end
end

"""
    getnode(root, s)

Return the node of the [`SimplexTree`](@ref) rooted at `root`
representing the simplex `s`.
If `s` is not in the complex,
return nothing.
"""
function getnode(complex::SimplexTree, simplex)
    if isempty(simplex)
        return complex
    elseif !haschild(complex, simplex[1])
        return nothing
    else
        return getnode(complex.children[simplex[1]], simplex[2:end])
    end
end

# The parent path from a node to the root represents a simplex
function simplex_from_node(node::SimplexTree)
    simplex = Array{Int, 1}()
    while !isemptyface(node)
        pushfirst!(simplex, node.label)
        node = node.parent
    end
    return simplex
end

## Modifying a complex

"""
    insert!(complex, simplex)
Insert a simplex and all its faces on an existing simplicial complex.
"""
@inline function insert!(complex::SimplexTree, simplex)
    return insert_ordered!(complex, sort(simplex))
end

function insert_ordered!(complex::SimplexTree, simplex)
    for (index, vertex) in enumerate(simplex)
        if !haschild(complex, vertex)
            push!(complex.children, STNode(vertex, complex))
        end
        insert_ordered!(child_from_label(complex,vertex), simplex[index+1:end])
    end
    return complex
end

function insert_ordered_face!(st::SimplexTree, simplex)
    node = st
    for vertex in simplex
        if haschild(node, vertex)
            node = child_from_label(node, vertex)
        else
            newnode = STNode(vertex, node)
            push!(node.children, newnode)
            node = newnode
        end
    end
    return st
end

## Folds

"""
    foldST(st::SimplexTree; basecase, joiner)

Fold st applying the function `basecase` to the leafs
and joining siblings using the function `joiner`.
"""
function foldST(node; basecase, joiner)
    if isleaf(node)
        return basecase(node)
    else
        return joiner(map(x -> foldST(x; basecase=basecase, joiner=joiner), node.children))
    end
end

"""
    foldST_depth(st::SimplexTree, depth; basecase, joiner)

Fold `st` until a given `depth`
applying the function `basecase` to the leafs
and joining siblings vi a the function `joiner`.

See also [`foldST`](@ref).
"""
function foldST_depth(node, depth; basecase, joiner)
    if depth == 0
        return basecase(node)
    else
        return joiner(map(x -> foldST_depth(x, depth-1;
                                            basecase=basecase,
                                            joiner=joiner)
                          , node.children)
                     )
    end
end

end
