module SparseTensors

export SparseTensor

using Base: size, similar, getindex, setindex!

struct SparseTensor{T<:Number, n} <: AbstractArray{T,n}
    data :: Dict{NTuple{n,Int}, T}
    dims :: NTuple{n,Int}
end

function SparseTensor(::Type{T}, dims::Int...) where {T}
    return SparseTensor(T, dims)
end

function SparseTensor(::Type{T}, dims::NTuple{n,Int}) where {T,n}
    return SparseTensor{T,n}(Dict{NTuple{n,Int}, T}(), dims)
end

## Interface to AbstractArray
Base.size(x::SparseTensor) = x.dims

Base.similar(::SparseTensor, ::Type{T}, dims::Dims) where {T} = SparseTensor(T, dims)

Base.getindex(x::SparseTensor{T,n}, I::Vararg{Int,n}) where {T,n} = get(x.data, I, zero(T))

function Base.setindex!(x::SparseTensor{T,n}, v, I::Vararg{Int, n}) where {T,n}
    if v == zero(T)
        delete!(x.data, I)
    else
        x.data[I] = v
    end
    return v
end

# Linear Algebraic operations
function Base.:+(x::SparseTensor{T,n}, y::SparseTensor{T,n}) where {T,n}
    return SparseTensor{T,n}(merge(+, x.data, y.data), x.dims)
end

function Base.:-(x::SparseTensor{T,n}, y::SparseTensor{T,n}) where {T,n}
    return SparseTensor{T,n}(merge(-, x.data, y.data), x.dims)
end

function Base.:*(a::T, x::SparseTensor{T,n}) where {T<:Number,n}
    return map(t -> a*t, x)
end

Base.:*(x::SparseTensor{T,n}, a::T) where {T<:Number,n} = a*x

Base.zero(x::SparseTensor{T,n}) where {T,n} = SparseTensor(T, x.dims)

end #module
