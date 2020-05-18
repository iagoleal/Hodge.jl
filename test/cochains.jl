using Hodge
using Test
using Random

@testset "Cochains" begin
    @testset "Assignment and basic unary ops" begin
        K = SimplicialComplex([[2,4,5], [1,2,3,4]])

        f = Cochain(Int, K, 1)
        # Test basic acessor functions
        @test basering(f) === Int
        @test degree(f) === 1
        @test basespace(f) === K
        # Forms over Int cannot accept Float64 values
        @test_throws InexactError f[3,4] = 10.5
        @test_throws InexactError f[3,4] = 7 // 2
        # Unless conversion is possible
        @test (f[3,4] = 10.0) === 10.0
        # Permutations flip sign
        @test f[3,4] === -f[4,3]
        # Also reassignment
        f[4,3] = 100
        @test f[3,4] == -100
        # Check whether inexistent edges return zero
        f[3,5] = 1000
        @test f[3,5] === 0
        # Assigning zero should not store anything in the internal Dict
        f[1,2] = 0
        @test !haskey(f.data, (1,2))
        f[1,2] = 5
        @test haskey(f.data, (1,2))
        f[1,2] = 0
        @test !haskey(f.data, (1,2))
        # Out of vertices:
        @test iszero(f[1,6])
        @test iszero(f[0,1])

        # Scaling / Self-arithmetic
        @test_nowarn 3 * f
        @test_nowarn 3.0 * f
        @test_throws InexactError 3.5 * f
        @test iszero(0 * f)
        @test 1 * f == f
        @test 5 * f == f * 5
        @test f + zero(f) == f
        @test zero(f) + f == f
        @test f - zero(f) == f
        @test f - f == zero(f)

        # Test copies / conversions
        @inferred collect(f)
        @inferred copy(f)
        g = similar(f, Float64, 5)
        @test basespace(f) == basespace(g)
    end
    @testset "Cochains interacting" begin
        K = SimplicialComplex([[2,4,5], [1,2,3,4]])
        f = Cochain(Int, K, 2)
        f[1,2,3] = 10
        f[3,1,4] = 5
        f[2,4,5] = 9

        # All binary operation must fail for different basespaces
        g1 = Cochain(Int, SimplicialComplex(), 2)
        @test basespace(f) != basespace(g1)
        @test_throws AssertionError f + g1
        @test_throws AssertionError f - g1
        @test_throws AssertionError g1 + f
        @test_throws AssertionError g1 - f
        @test_throws AssertionError cup(f, g1)

        # All binary operations should fail for different base rings
        g2 = similar(f, Float64, degree(f))
        @test basespace(f) === basespace(g2)
        @test basering(f) !== basering(g2)
        @test_throws MethodError f + g2
        @test_throws MethodError f - g2
        @test_throws MethodError g2 + f
        @test_throws MethodError g2 - f
        @test_throws MethodError cup(f, g2)

        # All binary operations should fail for different degrees, except for cup
        g3 = similar(f, basering(f), degree(f)+1)
        @test basespace(f) === basespace(g3)
        @test basering(f) === basering(g3)
        @test_throws MethodError f + g3
        @test_throws MethodError f - g3
        @test_throws MethodError g3 + f
        @test_throws MethodError g3 - f
        cup(f, g3)

        # For same basespace, basering and degree:
        g = copy(f)
        g[1,2,3] = 14245
        g[4,5,2] = -17
        g[1,2,4] = 2344
        f + g
        g + f
        f - g
        g - f
        cup(f,g)
    end
    @testset "Type conversions" begin
        K = SimplicialComplex(((1,2,3), (3,4), (4,5), (5,6,1), (1,2,3,4,5)))
        dim = rand(0:4)
        f = Cochain(Int, K, dim)
        for s in simplices(K,dim)
            f[s...] = rand(Int)
        end
        @test float(f) isa Cochain{Float64}
        @test complex(f) isa Cochain{Complex{Int}}
        @test complex(float(f)) isa Cochain{Complex{Float64}}

        dim = rand(0:4)
        g = Cochain(Float64, K, dim)
        for s in simplices(K,dim)
            g[s...] = rand(Float64)
        end
        @test rationalize(g) isa Cochain{Rational{Int}}
        @test rationalize(BigInt, g) isa Cochain{Rational{BigInt}}
    end
    @testset "Hodge Decomposition" begin
        K = SimplicialComplex(((1,2,3), (3,4), (4,5), (5,6,1)))
        f = Cochain(Float64, K, 1)
        for s in simplices(K,1)
            f[s...] = rand()
        end
        a, b, c = hodge(f, atol=1e-9)

        # Test if reconstruction goes well
        result = coboundary(a) + coboundary_adj(b) + c
        @test norm(f - result, 2) ≈ 0 atol=1e-5
        @test norm(f - result, Inf) ≈ 0 atol=1e-5

        # Usual Tests for decompositon
        @test norm(laplacian(c), 2) ≈ 0 atol=1e-7
    end
    @testset "Harmonic forms" begin
        # Hexagon missing an inner edge
        H = SimplicialComplex(((1,2,3), (3,4), (4,5), (5,6,1)))

        # Define a known harmonic form
        h = Cochain(Float64, H, 1)
        h[1,2] =  1
        h[1,3] =  2
        h[1,5] = -2
        h[1,6] = -1
        h[2,3] =  1
        h[3,4] =  3
        h[4,5] =  3
        h[5,6] =  1
        # Test if d, δ, and Δ of it are properly returning zero
        @test iszero(coboundary(h))
        @test iszero(coboundary_adj(h))
        @test iszero(laplacian(h))

        # Test if the hodge decomposition works
        a, b, c = hodge(h)
        @test norm(a, Inf)   ≈ 0 atol=1e-6
        @test norm(b, Inf)   ≈ 0 atol=1e-6
        @test norm(h-c, Inf) ≈ 0 atol=1e-6
    end
end
