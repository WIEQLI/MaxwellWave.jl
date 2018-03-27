export rqi!

function rqi!(A::AbsMat{T},  # system matrix
              μ::T,  # guess eigenvalue
              x::AbsVec{T},  # guess eigenvector; stores solution eigenvector
              maxit::Integer=10,  # max number of Rayleight quotient iteration
              τ::Real=Base.rtoldefault(T)  # tolerance in solution
             ) where {T<:Union{Float,CFloat}}
    I = similar(A)
    fill!(I, 0.0)
    I[diagind(I)] = 1.0

    B = similar(A)
    xold = similar(x)
    ps = PardisoSolver();

    λ = μ
    n = 0
    while n < maxit
        n += 1

        xold .= x

        # Solve a linear system.

        # Version 1
        # x .= (A - μ.*I) \ xold

        # Version 2
        B .= A
        B -= μ.*I
        x .= B \ xold

        # Version 3
        B .= A
        B -= μ.*I
    	# set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)
        # set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM);
        # set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON);
        # set_solver!(ps, Pardiso.DIRECT_SOLVER);
        # pardisoinit(ps);
        solve!(ps, x, B, xold);
        # info("Pardiso performed $(get_iparm(ps, 7)) iterative refinement steps");

        normalize!(x)

        # info("λ = $λ")
        λ = x⋅(A*x)
        if abs(1 - μ/λ) ≤ τ
            break
        end
        μ = λ
    end

    return λ
end
