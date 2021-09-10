using FdeSolver
using Test
using SpecialFunctions
using MittagLeffler

@testset "FdeSolver.jl" begin

    tSpan = [0, 1.5]
    y0 = 0
    β = 0.9

    function F(t, n, β, y)

        return (40320 ./ gamma(9 - β) .* t[n] .^ (8 - β) .- 3 .* gamma(5 + β / 2) ./ gamma(5 - β / 2) .* t[n] .^ (4 - β / 2) .+ 9/4 * gamma(β + 1) .+ (3/2 .* t[n] .^ (β / 2) .- t[n] .^ 4) .^ 3 .- y[n] .^ (3/2))

    end

    t, Yapp = FDEsolver(F, tSpan, y0, β)
    
    @test @isdefined(t)
    @test @isdefined(Yapp)

end
