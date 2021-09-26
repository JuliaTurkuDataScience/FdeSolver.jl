using FdeSolver
using Test
using SpecialFunctions

@testset "FdeSolver.jl" begin

    tSpan = [0, 5]
    y0 = [1, 0.5, 0.3]
    β = [0.5, 0.2, 0.6]

    function F(t, n, β, y)

        F1 = 1 / sqrt(pi) * (((y[n, 2] - 0.5) * (y[n, 3] - 0.3))^(1 / 6) + t[n]^(1 / 2))
        F2 = gamma(2.2) * (y[n, 1] - 1)
        F3 = gamma(2.8) / gamma(2.2) * (y[n, 2] - 0.5)

        return [F1, F2, F3]

    end

    t, Yapp = FDEsolver(F, tSpan, y0, β, nothing)

    @test @isdefined(t)
    @test @isdefined(Yapp)

end
