using FdeSolver
using Plots

## inputs
tSpan = [0, 10]     # [intial time, final time]
x0 = [1, 1]             # intial value
β = 1            # order of the derivative
par = [16.0, 4.0] # [spring constant for a mass on a spring, inertial mass]
h = 0.1

## Equation
function F(t, n, β, x, par)

      K = par[1]
      m = par[2]

      K ./ m .* x[n]

end

## Numerical solution
t, Yapp = FDEsolver(F, tSpan, x0, β, nothing, par, h = h)

#plot
plot(t, Yapp, linewidth = 5, title = "Solution of a 1D with order > 1",
     xaxis = "Time (t)", yaxis = "x(t)", label = "Approximation")
a = x0[1] .* map(cos, sqrt(par[1] / par[2]) .* t) .+ x0[2] ./ sqrt(par[1] / par[2]) .* map(sin, sqrt(par[1] / par[2]) .* t)
plot!(t, a, lw = 3, ls = :dash, label = "Exact solution")
