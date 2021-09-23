using FdeSolver
using Plots
using SpecialFunctions

## inputs
tSpan = [0, 1]     # [intial time, final time]
y0 = 0             # intial value
β = 0.9            # order of the derivative

# Equation
F(t, n, β, y) = (40320 ./ gamma(9 - β) .* t[n] .^ (8 - β) .- 3 .* gamma(5 + β / 2)
           ./ gamma(5 - β / 2) .* t[n] .^ (4 - β / 2) .+ 9/4 * gamma(β + 1) .+
           (3 / 2 .* t[n] .^ (β / 2) .- t[n] .^ 4) .^ 3 .- y[n] .^ (3 / 2))

# Jacobian
JacobF(t, n, β, y) = -(3 / 2) .* y[n] .^ (1 / 2)

## Numerical solution
t, Yapp = FDEsolver(F, tSpan, y0, β)

t1, Yapp1 = FDEsolver(F, tSpan, y0, β, JacobF)

#plot
plot(t, Yapp, linewidth = 5, title = "Solution of a 1D fractional IVP",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")
plot!(t1, Yapp1, linewidth = 5, ls = :dot, label = "Approximation with Jacobian")
plot!(t, t -> (t.^8 - 3 * t .^ (4 + β / 2) + 9 / 4 * t.^β),
      lw = 3, ls = :dash, label = "Exact solution")
