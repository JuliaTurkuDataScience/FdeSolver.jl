using Revise
push!(LOAD_PATH, "./examples")
using fdeSolver
using Plots

N = 1
I0 = 0.01

tSpan = [0, 100]
y0 = [1 - 2* I0, I0, I0]
α = .8*[1, 1, 1]
h = 0.1

function F(t, n, α, y)

    β = 0.4
    γ = 0.04

    dSdt = - β .* y[n, 1] .* y[n, 2]
    dIdt = β .* y[n, 1] .* y[n, 2] .- γ .* y[n, 2]
    dRdt = γ .* y[n, 2]

    return [dSdt, dIdt, dRdt]

end

t, Yapp = improveit8(F, tSpan, y0, α, h = h)

plot(t, Yapp, linewidth = 5, title = "Numerical solution of SIR model",
     xaxis = "Time (t)", yaxis = "y(t)")
