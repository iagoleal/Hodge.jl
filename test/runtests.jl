using Hodge
using Test

filepath(x) = joinpath(dirname(@__FILE__), x)

@testset "Hodge.jl" begin
    @info "Testing Simplicial Complexes"
    include(filepath("complexes.jl"))

    @info "Testing Cochains"
    include(filepath("cochains.jl"))
end
