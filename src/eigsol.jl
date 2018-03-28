export rqi!

function rqi!(A::AbsMat{T},  # system matrix
              μ::T,  # guess eigenvalue
              x::AbsVec{T},  # guess eigenvector; stores solution eigenvector
              maxit::Integer=10,  # max number of Rayleight quotient iteration
              τ::Real=Base.rtoldefault(T)  # tolerance in solution
             ) where {T<:Union{Float,CFloat}}
    m = size(A,1)
    I = similar(A)
    fill!(I, 0.0)
    I[diagind(I)] = 1.0

    B = similar(A)
    xold = similar(x)
    ps = PardisoSolver();
    set_msglvl!(ps, Pardiso.MESSAGE_LEVEL_ON)

    # First set the matrix type to handle general real symmetric matrices
    if T<:Real
        set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    else  # T<:Complex
        set_matrixtype!(ps, Pardiso.COMPLEX_NONSYM)
    end

    # Initialize the default settings with the current matrix type
    pardisoinit(ps)

    λ = μ
    n = 0
    while n < maxit
        n += 1

        xold .= x

        # Solve a linear system.

        # Version 1
        # x .= (A - μ.*I) \ xold

        # Version 2
        # B .= A
        # B -= μ.*I
        # x .= B \ xold

        # Version 3
        B .= A
        B -= μ.*I

        # Get the correct matrix to be sent into the pardiso function.
        # :N for normal matrix, :T for transpose, :C for conjugate
        A_pardiso = get_matrix(ps, B, :N)

        # Analyze the matrix and compute a symbolic factorization.
        set_phase!(ps, Pardiso.ANALYSIS)
        set_perm!(ps, randperm(m))
        pardiso(ps, A_pardiso, xold)

        # Compute the numeric factorization.
        set_phase!(ps, Pardiso.NUM_FACT)
        pardiso(ps, A_pardiso, xold)

        # Compute the solutions X using the symbolic factorization.
        set_phase!(ps, Pardiso.SOLVE_ITERATIVE_REFINE)
        # set_solver!(ps, Pardiso.ITERATIVE_SOLVER)
        pardiso(ps, x, A_pardiso, xold)

        normalize!(x)

        # info("λ = $λ")
        λ = x⋅(A*x)
        if abs(1 - μ/λ) ≤ τ
            break
        end
        μ = λ
    end

    # Free the PARDISO data structures.
    set_phase!(ps, Pardiso.RELEASE_ALL)
    pardiso(ps)

    return λ
end
