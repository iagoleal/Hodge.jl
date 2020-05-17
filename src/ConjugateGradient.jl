"""
    conjugate_gradient(A, b, inner; maxiters, atol, x0)

__Iteratively__ solves the linear equation `A(z) = b`,
where `A` is a linear operator and `b` an element of an abstract vector space.

This method does not assume that `A` is stored as a matrix
nor that `b` is stored as a vector,
only that the application `x -> A(x)` is well-defined
and that `A` represents a linear operator.

The optional parameters are:
    - `x0` is the initial value for the iterative method and defaults to zero.
    -  `atol` is the absolute tolerance of the stopping criterion.
        The method runs while the squared norm of the residue is larger than `atol`.
    - `maxiters` is the maximum number of iterations to do before stopping. Default is infinity.

Requirements:
    - `A` must be self-adjoint with respect to `inner`;
    - `A` must be positive-semidefinite.
    - If given, `x0` must be in the kernel of `A`, i.e. `A(x0) == zero(x0)`.

"""
function conjugate_gradient(A, b, inner; maxiters=nothing,
                                         atol=1e-6,
                                         x0=nothing)
    # Set number of iterations
    if isnothing(maxiters)
        iter = Iterators.cycle(1)
    else
        iter = 1:maxiters
    end
    # Set initial values
    if isnothing(x0)
        x = zero(b)
        r = b
    else
        x = x0
        r = b - A(x0)
    end
    p = r
    norm_r_prev = inner(r,r)
    for n in iter
        Ap    = A(p)
        alpha = norm_r_prev / inner(p, Ap)
        x     = x + alpha * p
        r     = r - alpha * Ap
        norm_r_new = inner(r,r)
        sqrt(norm_r_new) < atol && break

        beta  = norm_r_new / norm_r_prev
        p     = r + beta * p
        norm_r_prev = norm_r_new
    end
    return x
end
