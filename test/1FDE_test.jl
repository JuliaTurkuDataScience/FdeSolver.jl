using FdeSolver
using Test
using SpecialFunctions
using Statistics

@testset "FdeSolver.jl" begin

    tSpan = [0, 1]
    y0 = 0
    β = 0.9

    par = β
    F(t, y, par) = (40320 ./ gamma(9 - par) .* t .^ (8 - par) .- 3 .* gamma(5 + par / 2)
               ./ gamma(5 - par / 2) .* t .^ (4 - par / 2) .+ 9/4 * gamma(par + 1) .+
               (3 / 2 .* t .^ (par / 2) .- t .^ 4) .^ 3 .- y .^ (3 / 2))

    JacobF(t, y, β) = -(3 / 2) .* y .^ (1 / 2)

    t, Yapp = FDEsolver(F, tSpan, y0, β, par, StopIt = "Convergence", tol = 10e-8, itmax = 15)
    t1, Yapp1 = FDEsolver(F, tSpan, y0, β, par, JF = JacobF, StopIt = "Convergence", tol = 10e-8, itmax = 15)

    @test @isdefined(t)
    @test @isdefined(Yapp)

    @test @isdefined(t1)
    @test @isdefined(Yapp1)

    # @test abs(mean(Yapp) - 0.99971) < 0.05
    # @test abs(mean(Yapp1) - 0.99971) < 0.05

end
