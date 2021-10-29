using FdeSolver
using Test
using Statistics

@testset "FdeSolver.jl" begin

tSpan = [0, 10]     # [intial time, final time]
x0 = [1, 1]             # intial value
β = 2            # order of the derivative
par = [16.0, 4.0] # [spring constant for a mass on a spring, inertial mass]
h = 0.01

function F(t, x, par)

      K = par[1]
      m = par[2]

      - K ./ m .* x

end

t, Yapp = FDEsolver(F, tSpan, x0, β, par, h = h)

a = x0[1] .* map(cos, sqrt(par[1] / par[2]) .* t) .+ x0[2] ./ sqrt(par[1] / par[2]) .* map(sin, sqrt(par[1] / par[2]) .* t)

@test @isdefined(t)
@test @isdefined(Yapp)
@test abs((mean(Yapp) - mean(a))) < 0.05

end
