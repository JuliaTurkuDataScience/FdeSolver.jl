using FdeSolver
using Plots

## inputs
tSpan = [0, 25]        # [intial time, final time]
y0 = [34, 6]           # initial values
β = [0.98, 0.99]       # order of derivatives

## ODE model

par = [0.55, 0.028, 0.84, 0.026] # model parameters

function F(t, n, β, y, par)

    α1 = par[1]     # growth rate of the prey population
    β1 = par[2]     # rate of shrinkage relative to the product of the population sizes
    γ = par[3]      # shrinkage rate of the predator population
    δ = par[4]      # growth rate of the predator population as a factor of the product
                    # of the population sizes

    u = y[n, 1]      # population size of the prey species at time t[n]
    v = y[n, 2]      # population size of the predator species at time t[n]

    F1 = α1 .* u .- β1 .* u .* v
    F2 = - γ .* v .+ δ .* u .* v

    [F1, F2]

end

## Solution
t, Yapp = FDEsolver(F, tSpan, y0, β, par)

# plotting
plot(t, Yapp, linewidth = 5, title = "Solution to LV model with 2 FDEs",
     xaxis = "Time (t)", yaxis = "y(t)", label = ["Prey" "Predator"])
plot!(legendtitle="Population of")
