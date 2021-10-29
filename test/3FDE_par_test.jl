using FdeSolver
using Test
using SpecialFunctions

@testset "FdeSolver.jl" begin

    tSpan = [0, 5]
    y0 = [1, 0.5, 0.3]
    β = [0.5, 0.2, 0.6]

    function F(t, y)

        F1 = 1 / sqrt(pi) * (((y[ 2] - 0.5) * (y[3] - 0.3))^(1 / 6) + t^(1 / 2))
        F2 = gamma(2.2) * (y[1] - 1)
        F3 = gamma(2.8) / gamma(2.2) * (y[2] - 0.5)

        return [F1, F2, F3]

    end

    t, Yapp = FDEsolver(F, tSpan, y0, β)

    @test @isdefined(t)
    @test @isdefined(Yapp)

end
