using FdeSolver
using Plots
using SpecialFunctions

## inputs
tSpan = [0, 1]     # [intial time, final time]
y0 = 0             # intial value
β = 0.5            # order of the derivative
h=2^-5

# Equation
par = β
F(t, y, par) = (40320 ./ gamma(9 - par) .* t .^ (8 - par) .- 3 .* gamma(5 + par / 2)
           ./ gamma(5 - par / 2) .* t .^ (4 - par / 2) .+ 9/4 * gamma(par + 1) .+
           (3 / 2 .* t .^ (par / 2) .- t .^ 4) .^ 3 .- y .^ (3 / 2))
# Jacobian
JacobF(t, y, par) = -(3 / 2) .* y .^ (1 / 2)

## Numerical solution
t, Yapp = FDEsolver(F, tSpan, y0, β, par, h = h, StopIt = "Convergence", tol = 10e-8, itmax = 30)

t1, Yapp1 = FDEsolver(F, tSpan, y0, β, par, JF = JacobF, tol = 10e-8, itmax = 30)

#plot
plot(t, Yapp, linewidth = 5, title = "Solution of a 1D fractional IVP",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")
plot!(t1, Yapp1, linewidth = 5, ls = :dot, label = "Approximation with Jacobian")
plot!(t, t -> (t.^8 - 3 * t .^ (4 + β / 2) + 9 / 4 * t.^β),
      lw = 3, ls = :dash, label = "Exact solution")
