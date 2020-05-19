using Hodge
using Test
using Random

@testset "Simplicial Complexes" begin
    @testset "Construction" begin
        dim = 3
        simplex = rand(Int, dim+1)
        @test_throws UndefKeywordError B = SimplicialComplex(simplex)
        B = SimplicialComplex([simplex])
        for k in 0:dim
            Bk = skeleton(B, k)
            # Until k, has the same simplices
            for i in 0:k
                ss = simplices(Bk, i)
                @test simplices(B, i) == ss
                @test !isempty(ss)
            end
            # After that, is empty
            for i in k+1:dim
                ss = simplices(Bk, i)
                @test simplices(B, i) != ss
                @test isempty(ss)
            end
        end
    end
    @testset "Topological Operators" begin
        # Decomposition of sphere S^2
        S = skeleton(SimplicialComplex([(1,2,3,4)]), 2)
        @test euler_characteristic(S) == 2
        @test betti(S) == [1,0,1]

        # Decomposition of Torus T^2
        vs = [1 4 5 1; 3 8 9 3; 2 6 7 2; 1 4 5 1]
        tris = NTuple{3,Int}[]
        for (i,j) in Iterators.product(1:3, 1:3)
            push!(tris, (vs[i,j], vs[i,j+1], vs[i+1,j+1]))
            push!(tris, (vs[i,j], vs[i+1,j], vs[i+1,j+1]))
        end
        T = SimplicialComplex(tris)
        @test euler_characteristic(T) == 0
        @test betti(T) == [1,2,1]

        # Disjoint sum of two spheres S^2 ̧∐ S^2
        SS = SimplicialComplex([[2,3,4],
                                [1,3,4],
                                [1,2,4],
                                [1,2,3],
                                [6,7,8],
                                [5,7,8],
                                [5,6,8],
                                [5,6,7]
                               ])
        @test euler_characteristic(SS) == 4
        @test betti(SS) == [2,0,2]
    end
end
