import Combinatorics: levicivita, permutations
import SparseArrays
import LinearAlgebra

using Base: +, *, -, zero, iszero, ==
using Base: size, similar, getindex, setindex!
using Base: float, complex, rationalize
using Base: show
import Base.Iterators: take, drop

"""
    Cochain{R, n}

Represent the `n`-th group ``C^n(K; R)``
of cochains over the ring `R`
whose basespace is the simplicial complex ``K``.

For constructing Cochains,
see also the methods [`zero_cochain`](@ref)
and [`indicator_cochain`](@ref).

The elements of this type may be seem as
functions from the n-simplices of `K` to `R`
or as skew-symmetric n-tensors over the vertices of `K`.
This second perspective follows the ideas from the paper:
- Jiang, X., Lim, L., Yao, Y. et al. Statistical ranking and combinatorial Hodge theory. Math. Program. 127, 203–244 (2011). https://doi.org/10.1007/s10107-010-0419-x
"""
struct Cochain{R<:Number, n}
    basespace :: SimplicialComplex
    data      :: Dict{Tuple{Vararg{Int}}, R}
    function Cochain{R, n}(K::SimplicialComplex) where {R,n}
        return new(K, Dict{Tuple{Vararg{Int}}, R}())
    end
end

"""
    Cochain{R, n}(K::SimplicialComplex, itr)

Construct a Cochain whose basespace is `K`
and whose values are filled according to `pairs(itr)`.

#Example
```jldocstring
```
"""
function Cochain{R, n}(K::SimplicialComplex, itr) where {R,n}
    f = Cochain{R,n}(K)
    for (k,v) in pairs(itr)
        f[k...] = v
    end
    return f
end

"""
    zero_cochain(R, K, n)

Construct an identically zero `n`-cochain over `R`
and whose base space is `K`.
"""
@inline function zero_cochain(::Type{R}, K::SimplicialComplex, n::Int) where {R}
    return Cochain{R, n}(K)
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
    f = zero_cochain(R, K, simplex_dim(simplex))
    f[simplex...] = one(R)
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
    return Cochain{T, m}(basespace(f))
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

"""
    vectorfy(ω::Cochain)

Return an iterator over the values of a [`Cochain`](@ref)
viewing it as a vector.
"""
@inline vectorfy(f::Cochain) = (f[s...] for s in simplices(basespace(f), degree(f)))

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

@inline function Base.zero(f::Cochain)
    return zero_cochain(basering(f), basespace(f), degree(f))
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

Base.iszero(f::Cochain) = all(iszero, vectorfy(f))

function Base.:(==)(f::Cochain{R,n}, g::Cochain{R,n}) where {R,n}
    assert_basespaces(f, g)
    return iszero(f - g)
end

"""
    norm(ω[, p])

Calculate the p-norm of the [`Cochain`](@ref) `ω`.

By default, `p=2`.

Return a floating point,
no matter the base ring of `ω`.
"""
function norm(f::Cochain, p::Real=2) :: Float64
    if isinf(p)
        return maximum(abs, vectorfy(f))
    else
        return sum(t -> abs(t)^p, vectorfy(f)) ^ inv(p)
    end
end

"""
    norm2(ω)

Calculate the square of the usual inner product norm of a [`Cochain`](@ref) `ω`.
"""
norm2(f::Cochain) = sum(abs2, vectorfy(f))

@doc raw"""
    inner(ω, ξ)

Usual inner product between [`Cochain`](@ref)s.

!!! warning
    For complex Cochains,
    the conjugation is taken on the __first__ entry.

This inner product sees a n-cochain as a free vector space
over the (non-oriented) n-simplices of their base space.
Formally,
```math
\sum_{\sigma \in \mathrm{simplices}(K,n)} \overline{f(σ)} g(σ).
```
"""
function inner(f::Cochain{R,n}, g::Cochain{R,n}) :: R where {R,n}
    assert_basespaces(f,g)
    mutual_keys = intersect(keys(f.data), keys(g.data))
    if isempty(mutual_keys)
        return zero(R)
    else
        return sum(conj(g[i...]) * f[i...] for i in mutual_keys)
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
    h = zero_cochain(R, basespace(f), n + m)
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
    df = zero_cochain(basering(f), basespace(f), degree(f) + 1)
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
function coboundary_adj(f::Cochain)
    K  = basespace(f)
    δf = zero_cochain(basering(f), K, degree(f) - 1)
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
@inline function laplacian(f::Cochain)
    return coboundary(coboundary_adj(f)) + coboundary_adj(coboundary(f))
end

@doc raw"""
    hodge(ω)

Hodge decomposition of a [`Cochain`](@ref)
using [`inner`](@ref) as the inner product.

Return a tuple `(α,β,γ)` such that
```
ω == coboundary(α) + coboundary_adj(β) + γ
laplacian(γ) == 0
```
"""
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

#############################################
# Helper functions used throughout the code #
#############################################

# Helper function to calculate permutation sign
signperm = levicivita ∘ sortperm

# Interator yielding the (k-1)-faces of k-simplex
@inline faces(s) = ((s[j] for j in eachindex(s) if j != i) for i in eachindex(s))
@inline simplex_dim(l) = length(l) - 1

@inline function assert_basespaces(f,g)
    @assert(basespace(f) == basespace(g), "The Cochains are defined over different simplicial complexes")
end

function matrixify(L, K::SimplicialComplex, deg_from::Integer, deg_to::Integer=deg_from; sparse::Bool=true)
    n = numsimplices(K, deg_from)
    m = numsimplices(K, deg_to)
    A = sparse ? SparseArrays.spzeros(n, m) : zeros(n,m)
    for (i, s) in enumerate(simplices(K, deg_from))
        Ae_s = (vectorfy ∘ L ∘ indicator_cochain)(Float64, K, s)
        for (j,v) in enumerate(Ae_s)
            A[j, i] = v
        end
    end
    return A
end
