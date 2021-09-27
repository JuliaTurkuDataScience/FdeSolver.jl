using FdeSolver
using Test
using SpecialFunctions
using Statistics

@testset "FdeSolver.jl" begin

    tSpan = [0, 1.5]
    y0 = 0
    β = 0.9

    function F(t, n, β, y)

        return (40320 ./ gamma(9 - β) .* t[n] .^ (8 - β) .- 3 .* gamma(5 + β / 2) ./ gamma(5 - β / 2) .* t[n] .^ (4 - β / 2) .+ 9/4 * gamma(β + 1) .+ (3/2 .* t[n] .^ (β / 2) .- t[n] .^ 4) .^ 3 .- y[n] .^ (3/2))

    end

    JacobF(t, n, β, y) = -(3 / 2) .* y[n] .^ (1 / 2)

    t, Yapp = FDEsolver(F, tSpan, y0, β, nothing, StopIt = "Convergence", tol = 10e-8, itmax = 15)
    t1, Yapp1 = FDEsolver(F, tSpan, y0, β, JacobF, StopIt = "Convergence", tol = 10e-8, itmax = 15)

    @test @isdefined(t)
    @test @isdefined(Yapp)

    @test @isdefined(t1)
    @test @isdefined(Yapp1)

    @test abs(mean(Yapp) - 0.99971) < 0.05
    @test abs(mean(Yapp1) - 0.99971) < 0.05

end
