using Test
using Hodge

# Define an known harmonic form
# and test if d, δ, and Δ of it are properly returning zero
H = SimplicialComplex()
insert!(H, (1,2,3))
insert!(H, (3,4))
insert!(H, (4,5))
insert!(H, (5,6,1))

h = Cochain(Float64, H, 1)
h[1,2] =  1
h[1,3] =  2
h[1,5] = -2
h[1,6] = -1
h[2,3] =  1
h[3,4] =  3
h[4,5] =  3
h[5,6] =  1

Δh = laplacian(h)
dh = coboundary(h)
δh = coboundary_adj(h)

# Test whether the operators return zero for an harmonic form
@test iszero(Δh)
@test iszero(dh)
@test iszero(δh)

# Test the hodge decomposition for harmonic forms
α, β, γ = hodge(h)
@test iszero(α)
@test iszero(β)
@test γ == h
