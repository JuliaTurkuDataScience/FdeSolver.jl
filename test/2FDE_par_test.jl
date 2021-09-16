using FdeSolver
using Test
using Statistics

@testset "FdeSolver.jl" begin

    tSpan = [0, 25]
    y0 = [34, 6]
    β = [0.98, 0.99]
    par = [0.55, 0.028, 0.80, 0.024]

    function F(t, n, β, y, par)

        F1 = par[1] .* y[n, 1] .- par[2] .* y[n, 1] .* y[n, 2]
        F2 = - par[3] .* y[n, 2] .+ par[4] .* y[n, 1] .* y[n, 2]

        [F1, F2]

    end

    t, Yapp = FDEsolver(F, tSpan, y0, β, par)

    @test @isdefined(t)
    @test @isdefined(Yapp)
    @test mean(t) == 12.5
    @test abs((mean(Yapp) - 28.192035534885537)) < 0.05

end
