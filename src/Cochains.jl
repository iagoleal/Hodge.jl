import Combinatorics: levicivita, permutations
using ComputedFieldTypes

using Base: +, *, -, zero, iszero, ==
using Base: size, similar, getindex, setindex!
using Base: show
import Base.Iterators: take, drop

# Helper function to calculate permutation sign
signperm = levicivita ∘ sortperm

# Interator yielding the (k-1)-faces of k-simplex
@inline faces(s) = ((s[j] for j in eachindex(s) if j != i) for i in eachindex(s))
@inline simplex_dim(l) = length(l) - 1

@inline function assert_basespaces(f,g)
    @assert(basespace(f) == basespace(g), "The Cochains are defined over different simplicial complexes")
end

"""
    Cochain{R, n}

Represent the `n`-th group ``C^n(K; R)``
of cochains over the ring `R`
whose basespace is the simplicial complex ``K``.

The elements of this type may be seem as
functions from the n-simplices of `K` to `R`
or as skew-symmetric n-tensors over the vertices of `K`.
This second perspective follows the ideas from the paper:
- Jiang, X., Lim, L., Yao, Y. et al. Statistical ranking and combinatorial Hodge theory. Math. Program. 127, 203–244 (2011). https://doi.org/10.1007/s10107-010-0419-x
"""
struct Cochain{R<:Number, n}
    basespace :: SimplicialComplex
    values    :: Dict{Tuple{Vararg{Int}}, R}
end

## Constructing Cochains

# Zero cochain over K
"""
    Cochain(R, K, n)

Construct an identically zero n-cochain over `R`
and whose base space is `K`.
"""
@inline function Cochain(::Type{R}, K::SimplicialComplex, n::Int) where {R}
    return Cochain{R, n}(K, Dict{Tuple{Vararg{Int}}, R}())
end

# Basis cochains
"""
    indicator_cochain(R, K, σ)

Return the indicator function `f` of the n-simplex `σ` as a [`Cochain`](@ref).
That is, a cochain such that `f(σ) = 1` and  `f(τ) = 0`
for all other n-simplices of `K`.

Notice that, nevertheless,
`f` is still skew-symmetric over permutations of `σ`s indices.

If `Κ` does not contain `σ`,
the returned cochain is identically zero.
"""
function indicator_cochain(::Type{R}, K::SimplicialComplex, simplex) where {R}
    f = Cochain(R, K, simplex_dim(simplex))
    if hassimplex(K, simplex)
        f[simplex...] = one(R)
    end
    return f
end

## Accessor functions
"""
    basespace(ω::Cochain)

Return the simplicial complex
that `ω` is associated with.
"""
@inline basespace(f::Cochain) = f.basespace

"""
    basering(ω::Cochain)

Return the Ring that `ω` is defined over.
"""

@inline basering(::Cochain{R,n}) where {R,n} = R
"""
    degree(ω::Cochain)

Return the degree of a [`Cochain`](@ref).
"""
@inline degree(::Cochain{R,n}) where {R,n} = n


## Interface for AbstractArray type
Base.IndexStyle(::Type{<:Cochain}) = IndexCartesian()

function Base.size(f::Cochain)
    return ntuple(i -> numvertices(basespace(f)), degree(f)+1)
end

function Base.similar(f::Cochain, ::Type{T}, m::Int) where {T}
    return Cochain(T, basespace(f), m)
end

function Base.getindex(f::Cochain{R, n}, I::Vararg{Int}) where {R,n}
    length(I) == degree(f) + 1 || throw(BoundsError(f, I))
    perm    = sortperm(collect(I))
    simplex = I[perm]
    return levicivita(perm) * get(f.values, simplex, zero(R))
end

function Base.setindex!(f::Cochain{R,n}, v, I::Vararg{Int}) where {R,n}
    length(I) == degree(f) + 1 || throw(BoundsError(f, I))
    perm    = sortperm(collect(I))
    simplex = I[perm]
    if hassimplex(basespace(f), simplex)
        if iszero(v)
            delete!(f.values, simplex)
        else
            f.values[simplex] = levicivita(perm) * v
        end
    end
    return v
end

function Base.copy(f::Cochain{R,n}) where {R,n}
    return Cochain{R,n}(basespace(f), copy(f.values))
end

function Base.collect(f::Cochain{R,n}) where {R,n}
    A = zeros(R, size(f))
    for (k,v) in f.values
        for p in permutations(k)
            A[p...] = signperm(p) * v
        end
    end
    return A
end

## Showing Cochains on REPL
function Base.show(io::IO, f::Cochain)
    show(io, collect(f))
end

function Base.show(io::IO, ::MIME"text/plain", f::Cochain)
    #= show(io, "text/plain", collect(f)) =#
    println(io, degree(f), "-th degree Cochain over ", basering(f), ":")
    Base.print_array(io, collect(f))
end

## Vector Space structure

@inline function Base.:+(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f, g)
    return Cochain{R,n}(basespace(f), merge(+, f.values, g.values))
end

@inline function Base.:-(f::Cochain{R,n}) where {R,n}
    return Cochain{R,n}(basespace(f), Dict(k => -v for(k,v) in f.values))
end

@inline function Base.:-(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f, g)
    return f + (-g)
end

@inline function Base.zero(f::Cochain{R,n}) where {R,n}
    Cochain(R, basespace(f), n)
end

@inline function Base.:*(a::R, f::Cochain{R,n}) where {R<:Number,n}
    if iszero(a)
        return zero(f)
    else
        return Cochain{R, n}(basespace(f), Dict(k => a*v for (k,v) in f.values))
    end
end

@inline Base.:*(f::Cochain, a) = a * f

Base.iszero(f::Cochain) = all(x -> iszero(x.second), f.values)

function Base.:(==)(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f, g)
    return iszero(f - g)
end
"""
    norm(ω[, p])

Calculate the p-norm of the [`Cochain`](@ref) `ω`.

By default, `p=2`.
"""
norm(f::Cochain, p::Real=2) =
    isinf(p) ? maximum(collect(f)) : sum((t -> abs(t)^p).(f)) ^ inv(p)

"""
    norm2(ω)

Calculate the square of the usual inner product norm of a [`Cochain`](@ref) `ω`.
"""
norm2(f::Cochain) = sum(abs2.(f))

"""
    inner(ω, ξ)

Usual inner product between [`Cochain`](@ref)s.
"""
function inner(f::Cochain, g::Cochain)
    assert_basespaces(f,g)
    return sum(f .* g)
end

## Topological Operators

"""
    cup(ω, ξ)

The cup product ``\\omega \\smile \\xi``
of two [`Cochain`](@ref)s.
"""
function cup(f::Cochain{R,n}, g::Cochain{R,m}) where {R,n,m}
    assert_basespaces(f, g)
    h = Cochain(R, basespace(f), n + m)
    for s in simplices(basespace(f), n + m)
        h[s...] = f[take(s,n+1)...] * g[drop(s,n)...]
    end
    return h
end

"""
    coboundary(ω)

The coboundary or __discrete exterior derivative__
of a [`Cochain`](@ref).

The coboundary of ``ω`` applied to a simplex ``σ``
equals the alternating sum of ``ω`` applied to the faces of ``σ``.
"""
function coboundary(f::Cochain{R,n}) where {R,n}
    df = Cochain(R, basespace(f), n+1)
    for s in simplices(basespace(f), n+1)
        sg = 1
        for t in faces(s)
            df[s...] += sg * f[t...]
            sg = -sg
        end
    end
    return df
end

## Dependent on a choice of inner product

"""
    coboundary_adj(ω[, inner])

    The adjoint of the [`coboundary`](@ref)
with respect to a given inner product.

If the variable inner is omitted,
the default is the canonical inner product
treating `ω` as a n-th order tensor.
"""
function coboundary_adj(f::Cochain{R,n}, inner=inner) where {R,n}
    K = basespace(f)
    δf = Cochain(R, K, degree(f)-1)
    for s in simplices(K, degree(f)-1)
        e_s = (1/n) * indicator_cochain(R, K, s)
        δf[s...] = inner(f, coboundary(e_s))
    end
    return δf
end

"""
    laplacian(ω[, inner])

The (higher-order) laplacian of `ω`,
defined as
```math
    \\delta = d d^* + d^* d
```
where ``d`` and ``d^*`` are, respectively,
the [`coboundary`](@ref) and [`coboundary_adj`](@ref) operators.

The parameter `inner` dictates the inner product
needed for the laplacian's definition.
If it is omitted, the default is the usual dot product
treating `ω` as a n-th order tensor.
"""
@inline function laplacian(f::Cochain, inner=inner)
    return coboundary(coboundary_adj(f, inner)) + coboundary_adj(coboundary(f), inner)
end

@doc raw"""
    hodge(ω[, inner])

Hodge decomposition of a [`Cochain`](@ref)
using `inner` as the inner product.

Return a tuple `(α,β,γ)` such that
```
ω == coboundary(α) + coboundary_adj(β) + γ
laplacian(γ) == 0
```

The parameter `inner` dictates the inner product
needed for the `coboundary_adj` definition.
If it is omitted, the default is the usual dot product
treating `ω` as a n-th order tensor.
"""
function hodge(f::Cochain, inner=inner)
    L = x -> laplacian(x,inner) + x
    u = conjugate_gradient(L, f, inner)
    alpha = coboundary_adj(u, inner)
    beta  = coboundary(u)
    gamma = u
    return alpha, beta, gamma
end
#= The method used here is the following:
   We want to solve the following equation:
   ω = dα + dβ + γ,   where Δγ = 0.

   Lemma 1: C^n = im(Δ) ⊕ ker(Δ)
     Proof: Let x in im(Δ)  (i.e. x = Δu, for some u in C^k)
                y in ker(Δ) (i.e. Δy = 0)
            Then <x,y> = <Au,y> = <u,Ay> = <u,0> = 0.

   Therefore, there is γ harmonic such that
         ω = Δu + γ
         ω = d(δu) + δ(du) + γ.

   To solve this equation in a basis free manner,
   we use the conjugate gradient method.
   Since it requires the operator to be positive-definite,
   we supply L = Δ + λI to the method. (where λ > 0)
   This gives us a unique solution u satisfying
         (Δ + λI)u = ω
         ω = Δu + λu
         ω = d(δu) + δ(du) + λu.
   Notice that, by Lemma 1, λu must be harmonic.

   By setting:
         α = δu
         β = du
         γ = λu,
   We have three Cochains satisfying
         ω = dα + δβ + γ.
=#
