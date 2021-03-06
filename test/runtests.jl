using Hodge
using Test
using Random

filepath(x) = joinpath(dirname(@__FILE__), x)

# Fix random number seed, for reproducibility

@testset "Hodge.jl" begin
    @info "Testing Simplicial Complexes"
    include(filepath("complexes.jl"))

    @info "Testing Cochains"
    include(filepath("cochains.jl"))
end
