using FdeSolver
using Test
using Statistics

@testset "FdeSolver.jl" begin

    tSpan = [0, 25]
    y0 = [34, 6]
    β = [0.98, 0.99]
    par = [0.55, 0.028, 0.80, 0.024]

    function F(t, y, par)

        F1 = par[1] .* y[1] .- par[2] .* y[1] .* y[2]
        F2 = - par[3] .* y[2] .+ par[4] .* y[1] .* y[2]

        [F1, F2]

    end

    t, Yapp = FDEsolver(F, tSpan, y0, β, par)

    @test @isdefined(t)
    @test @isdefined(Yapp)
    @test abs((mean(Yapp) - 28.19203)) < 0.05

end
