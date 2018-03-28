@testset "eigsol" begin

function poisson2d(∆x, ∆y, Nx, Ny)
    N = Nx*Ny
    ∆x⁻², ∆y⁻² = 1.0/∆x^2, 1.0/∆y^2
    row₀ = reshape(1:N, (Nx,Ny)); col₀ = row₀; val₀ = 2.0(∆x⁻² + ∆y⁻²)

    rowx₊ = row₀[1:end-1,:]; colx₊ = circshift(col₀, (-1,0))[1:end-1,:]; valx₊ = -∆x⁻²
    rowx₋ = row₀[2:end,:]; colx₋ = circshift(col₀, (1,0))[2:end,:]; valx₋ = -∆x⁻²
    rowy₊ = row₀[:,1:end-1]; coly₊ = circshift(col₀, (0,-1))[:,1:end-1]; valy₊ = -∆y⁻²
    rowy₋ = row₀[:,2:end]; coly₋ = circshift(col₀, (0,1))[:,2:end]; valy₋ = -∆y⁻²

    A = sparse(row₀[:], col₀[:], val₀, N, N) +
        sparse(rowx₊[:], colx₊[:], valx₊, N, N) +
        sparse(rowx₋[:], colx₋[:], valx₋, N, N) +
        sparse(rowy₊[:], coly₊[:], valy₊, N, N) +
        sparse(rowy₋[:], coly₋[:], valy₋, N, N)

    return A
end

# Choose the target mode numbers
nx, ny = 2, 2  # mode number in x and y

# Construct the problem.
Nx, Ny = 99, 99
N = Nx*Ny
∆x = ∆y = 0.01
Lx, Ly = (Nx+1)*∆x, (Ny+1)*∆y
A = poisson2d(∆x, ∆y, Nx, Ny)

# Calculate the theoretical estimate of the eigenvalue.
λtheory = (π*nx/Lx)^2 + (π*ny/Ly)^2

# Perform Rayleigh quotient iteration to calculate the numerical eigenvalue.
x = rand(N)
tic()
λ = rqi!(A, λtheory, x)
toc()

# Prepare plot data.
u = zeros(Nx+2, Ny+2)
u[2:Nx+1, 2:Ny+1] .= reshape(x, (Nx,Ny))

xs = (0:Nx+1)*∆x .+ zeros(Ny+2)'
ys = (0:Ny+1)'*∆y .+ zeros(Nx+2)

λnumerical = λ
err = abs((λnumerical - λtheory)/λtheory)
# println("|λnumerical − λtheory| / |λtheory| = $err")
@test err ≤ 1e-3

end  # @testset "eigsol"
