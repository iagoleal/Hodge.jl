import .SparseTensors
import Combinatorics: levicivita

using Base: +, *, -, zero
using Base: size, similar, getindex, setindex!

# Helper function to calculate permutation sign
permsign = levicivita ∘ sortperm


@inline faces(s) = ((s[j] for j in eachindex(s) if j != i) for i in eachindex(s))

"""
    Cochain{R, K, n}

Represent the group ``C^n(K; R)``
of cochains over a simplicial complex ``K``
of dimension `n` and over the ring `R`.

The elements of this type may be seem as
functions from the n-simplices of `K` to `R`
or as skew-symmetric n-tensors over the vertices of `K`.
This second perspective follows the ideas from the paper:
- Jiang, X., Lim, L., Yao, Y. et al. Statistical ranking and combinatorial Hodge theory. Math. Program. 127, 203–244 (2011). https://doi.org/10.1007/s10107-010-0419-x
"""
struct Cochain{R, K, n} <: AbstractArray{R,n} where {R <: Number}
    values :: SparseTensor{R,n}
end

## Constructing Cochains

# Zero cochain over K
function Cochain(::Type{R}, K::SimplicialComplex, n::Int) where {R}
    bounds = ntuple(i -> numvertices(sc), n)
    return Cochain{R,K,n}(SparseTensor(R, bounds))
end

## Accessor functions
"""
    basespace(ω)

Return the simplicial complex
that `ω` is associated with.
"""
basespace(::Cochain{R,K,n}) where {R,K,n} = K

"""
    order(ω)
Return the order of a [`Cochain`](@ref)
"""
order(::Cochain{K,R,n}) where {R,K,n} = n

## Interface for AbstractArray type
Base.IndexStyle(::Type{<:Cochain}) = IndexCartesian()

Base.size(f::Cochain) = f.data.dims

Base.similar(::Cochain{R,K,n}, ::Type{T}, n) where {T,R,K,n}  = Cochain(T,K,n)

function Base.getindex(f::Cochain{R,K,n}, I::Vararg{Int,n}) where {R,K,n}
    order = sortperm(collect(I))
    return levicivita(order) * f.data[I[order]]
end

function Base.setindex!(f::Cochain{R,K,n}, v, I::Vararg{Int, n}) where {R,K,n}
    simplex = I[sortperm(collect(I))]
    if hassimplex(K, I)
        f.data[simplex] = v
    end
    return v
end

## Vector Space structure
function Base.:+(f::Cochain{R,K,n}, g::Cochain{R,K,n}) where {R,K,n}
    return Cochain{R,K,n}(f.data + g.data)
end

function Base.:-(f::Cochain{R,K,n}, g::Cochain{R,K,n}) where {R,K,n}
    return Cochain{R,K,n}(f.data - g.data)
end

function Base.:*(a::R, f::Cochain{R,K,n}) where {R,K,n}
    return Cochain{R,K,n}(a*f.data)
end

Base.:*(f::Cochain, a) = a*f

Base.zero(::Cochain{R,K,n})       where {R,K,n} = Cochain(R, K, n)
Base.zero(::Type{Cochain{R,K,n}}) where {R,K,n} = Cochain(R, K, n)

"""
    norm(ω[, p])

Calculate the p-norm of the [`Cochain`](@ref) `ω`.

The default is `p=2`.
"""
norm(f::Cochain, p=2) = sum(map(t -> abs(t)^p, f)) ^ inv(p)

"""
    norm2(ω)

Calculate the square of the usual inner product norm of a [`Cochain`](@ref) `ω`.
"""
norm2(f::Cochain) = sum(abs2, f)

"""
    inner(ω, ξ)

Usual inner product between [`Cochain`](@ref)s.
"""
inner(f::Cochain, g::Cochain) = sum(map(*, f, g))

## Topological Operators

@doc raw"""
    coboundary(ω)

The coboundary or __discrete exterior derivative__
of a [`Cochain`](@ref).

The coboundary of ``ω`` applied to a simplex ``σ``
equals the alternating sum of ``ω`` applied to the faces of ``σ``.
"""
function coboundary(f::Cochain{R,K,n}) where {R,K,n}
    df = Cochain{R,K, n+1}
    for s in simplices{k+1}
        sg = 1
        for t in faces(s)
            df[s...] += sg * f[t...]
            sg = -sg
        end
    end
    return df
end

@doc raw"""
    coboundary_adj(ω)

The adjoint of the coboundary
with respect to the usual inner product over the k-simplices
of a simplicial complex `K`.
This function satisfies
```math
    \langle ξ, \mathrm{coboundary_adj}(ω) \rangle_K
    = \langle boundary(ξ), ω \rangle_K
```
"""
function coboundary_adj(f::Cochain)
end

"""
    laplacian(ω)

The (higher-order) laplacian of `ω`,
defined as
```math
    \\delta = d d^* + d^* d
```
where ``d`` and ``d^*`` are, respectively,
the `coboundary` and `coboundary_adj` operators.
"""
@inline function laplacian(f::Cochain)
    return coboundary(coboundary_adj(f)) + coboundary_adj(coboundary(f))
end
