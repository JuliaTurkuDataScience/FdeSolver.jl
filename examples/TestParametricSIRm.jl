using Revise
push!(LOAD_PATH, "./examples")
using fdeSolver
using Plots

N = 1
I0 = 0.001

tSpan = [0, 100]
y0 = [1 - I0, I0, 0]
α = [1, 1, 1]
h = 0.1
par = [0.4, 0.04]

function F(t, n, α, y, par)

    dSdt = - par[1] .* y[n, 1] .* y[n, 2]
    dIdt = par[1] .* y[n, 1] .* y[n, 2] .- par[2] .* y[n, 2]
    dRdt = par[2] .* y[n, 2]

    return [dSdt, dIdt, dRdt]

end

t, Yapp = improveit10(F, tSpan, y0, α, par, h = h)

plot(t, Yapp, linewidth = 5, title = "Numerical solution of SIR model",
     xaxis = "Time (t)", yaxis = "y(t)")
