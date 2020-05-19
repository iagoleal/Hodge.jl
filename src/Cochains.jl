import Combinatorics: levicivita, permutations
import SparseArrays

using Base: +, *, -, zero, iszero, ==
using Base: size, similar, getindex, setindex!
using Base: float, complex, rationalize
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
    data      :: Dict{Tuple{Vararg{Int}}, R}
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


######################
# Accessor functions #
######################

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
    return levicivita(perm) * get(f.data, simplex, zero(R))
end

function Base.setindex!(f::Cochain{R,n}, v, I::Vararg{Int}) where {R,n}
    length(I) == degree(f) + 1 || throw(BoundsError(f, I))
    perm    = sortperm(collect(I))
    simplex = I[perm]
    if hassimplex(basespace(f), simplex)
        if iszero(v)
            delete!(f.data, simplex)
        else
            f.data[simplex] = levicivita(perm) * v
        end
    end
    return v
end

function Base.copy(f::Cochain{R,n}) where {R,n}
    return Cochain{R,n}(basespace(f), copy(f.data))
end

function Base.collect(f::Cochain{R,n}) where {R,n}
    A = zeros(R, size(f)) :: Array{R, n+1}
    for (k,v) in f.data
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

#########################
# Base Ring Conversions #
#########################

function Base.float(f::Cochain{R}) where {R}
    return Cochain{float(R), degree(f)}(basespace(f), Dict(k => float(v) for (k,v) in f.data))
end

function Base.complex(f::Cochain{R}) where {R}
    return Cochain{complex(R), degree(f)}(basespace(f), Dict(k => complex(v) for (k,v) in f.data))
end

function Base.rationalize(::Type{T}, f::Cochain{R}; kvs...) where {T<:Integer, R<:AbstractFloat}
    return Cochain{Rational{T}, degree(f)}(basespace(f), Dict(k => rationalize(T, v; kvs...) for (k,v) in f.data))
end

function Base.rationalize(f::Cochain{R}; kvs...) where {R<:AbstractFloat}
    return rationalize(Int, f; kvs...)
end

##########################
# Vector Space Structure #
##########################

@inline function Base.:+(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f, g)
    return Cochain{R,n}(basespace(f), merge(+, f.data, g.data))
end

@inline function Base.:-(f::Cochain{R,n}) where {R,n}
    return Cochain{R,n}(basespace(f), Dict(k => -v for(k,v) in f.data))
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
    elseif isone(a)
        return f
    else
        return Cochain{R, n}(basespace(f), Dict(k => a*v for (k,v) in f.data))
    end
end

function Base.:(*)(a::Number, f::Cochain{R}) where {R}
    return convert(R, a) * f
end

@inline Base.:*(f::Cochain, a) = a * f

Base.iszero(f::Cochain) = all(iszero, values(f.data))

function Base.:(==)(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f, g)
    return iszero(f - g)
end

"""
    norm(ω[, p])

Calculate the p-norm of the [`Cochain`](@ref) `ω`.

By default, `p=2`.
"""
function norm(f::Cochain, p::Real=2)
    if iszero(f)
        return zero(basering(f))
    elseif isinf(p)
        return maximum(values(f.data))
    else
        return sum(t -> abs(t)^p, values(f.data)) ^ inv(p)
    end
end

"""
    norm2(ω)

Calculate the square of the usual inner product norm of a [`Cochain`](@ref) `ω`.
"""
norm2(f::Cochain{R}) where R = iszero(f) ? zero(R) : sum(abs2, values(f.data))

@doc raw"""
    inner(ω, ξ)

Usual inner product between [`Cochain`](@ref)s.

This inner product sees a n-cochain as a free vector space
over the (non-oriented) n-simplices of their base space.
Formally,
```math
    \sum_{\sigma \in \mathrm{simplices}(K,n)} f(σ) g(σ).
```
"""
function inner(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f,g)
    mutual_keys = intersect(keys(f.data), keys(g.data))
    if isempty(mutual_keys)
        return zero(R)
    else
        return sum(g[i...] * f[i...] for i in mutual_keys)
    end
end


#########################
# Topological Operators #
#########################

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
function coboundary(f::Cochain)
    df = Cochain(basering(f), basespace(f), degree(f) + 1)
    for s in simplices(basespace(f), degree(f) + 1)
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
    coboundary_adj(ω)

The adjoint of the [`coboundary`](@ref)
with respect to usual inner product.
"""
#= If the variable inner is omitted, =#
#= the default is the canonical inner product =#
#= treating `ω` as a n-th order tensor. =#
function coboundary_adj(f::Cochain)
    K  = basespace(f)
    δf = Cochain(basering(f), K, degree(f) - 1)
    for s in simplices(K, degree(f) - 1)
        e_s = indicator_cochain(basering(f), K, s)
        δf[s...] = inner(f, coboundary(e_s))
    end
    return δf
end # Note: This only works for diagonal inner products...

"""
    laplacian(ω)

The (higher-order) laplacian of `ω`,
defined as
```math
    \\delta = d d^* + d^* d
```
where ``d`` and ``d^*`` are, respectively,
the [`coboundary`](@ref) and [`coboundary_adj`](@ref) operators.
"""
#= The parameter `inner` dictates the inner product =#
#= needed for the laplacian's definition. =#
#= If it is omitted, the default is the usual dot product =#
#= treating `ω` as a n-th order tensor. =#
@inline function laplacian(f::Cochain)
    return coboundary(coboundary_adj(f)) + coboundary_adj(coboundary(f))
end

@doc raw"""
    hodge(ω)

Hodge decomposition of a [`Cochain`](@ref)
using `inner` as the inner product.

Return a tuple `(α,β,γ)` such that
```
ω == coboundary(α) + coboundary_adj(β) + γ
laplacian(γ) == 0
```
"""
#= The parameter `inner` dictates the inner product =#
#= needed for the `coboundary_adj` definition. =#
#= If it is omitted, the default is the usual dot product =#
#= treating `ω` as a n-th order tensor. =#
function hodge(f::Cochain; atol=1e-9)
    maxiters = numsimplices(basespace(f), degree(f))
    λ = atol
    L = x -> laplacian(x) + λ * x
    u = conjugate_gradient(L, f, inner; atol=atol, maxiters=maxiters)
    alpha = coboundary_adj(u)
    beta  = coboundary(u)
    gamma = λ * u
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
   we supply L = Δ + εI to the method. (where ε > 0 is small)
   This gives us a unique solution u satisfying
         (Δ + εI)u = ω
         ω = Δu + εu
         ω = d(δu) + δ(du) + εu.
   By Lemma 1 and continuity,
   If ε is small enough, Δ(εu) ≈ 0.

   By setting:
         α = δu
         β = du
         γ = εu,
   We have three Cochains satisfying
         ω = dα + δβ + γ.
=#

##########################################
# Write cochains / operators as matrices #
##########################################

function vectorfy(f::Cochain{R, n}) where {R, n}
    return (f[s...] for s in simplices(basespace(f), n))
end

function matrixify(L, K::SimplicialComplex, deg_from::Integer, deg_to::Integer=deg_from)
    n = numsimplices(K, deg_from)
    m = numsimplices(K, deg_to)
    A = SparseArrays.spzeros(n, m)
    for (i, s) in enumerate(simplices(K, deg_from))
        Ae_s = (vectorfy ∘ L ∘ indicator_cochain)(Float64, K, s)
        for (j,v) in enumerate(Ae_s)
            A[j, i] = v
        end
    end
    return A
end
