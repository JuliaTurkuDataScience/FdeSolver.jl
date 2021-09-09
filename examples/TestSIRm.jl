using Revise
push!(LOAD_PATH, "./src")
using FdeSolver
using Plots

N = 1
I0 = 0.001

tSpan = [0, 100]
y0 = [1 - I0, I0, 0]
α = [1, 1, 1]
h = 0.1

function F(t, n, α, y)

    β = 0.4
    γ = 0.04

    dSdt = - β .* y[n, 1] .* y[n, 2]
    dIdt = β .* y[n, 1] .* y[n, 2] .- γ .* y[n, 2]
    dRdt = γ .* y[n, 2]

    return [dSdt, dIdt, dRdt]

end

t, Yapp = FDEsolver(F, tSpan, y0, α, h = h)

plot(t, Yapp, linewidth = 5, title = "Numerical solution of SIR model",
     xaxis = "Time (t)", yaxis = "y(t)")
