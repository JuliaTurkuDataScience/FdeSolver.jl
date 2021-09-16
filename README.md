# FdeSolver

[![CI](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl/branch/main/graph/badge.svg?token=SJ5F6RQ31P)](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl)

This is a Pkg in **Julia** for solution to a class of fractional differential equations and system equations.
Many advanced source codes are available in [MATLAB](https://www.dm.uniba.it/members/garrappa/software), but they are not open source projects like this one in Julia. Hence, the purpose is to develop a Julia package that numerically solves nonlinear fractional ordinary differential equations.

### Method

We implement the [predictor-corrector](https://link.springer.com/article/10.1023/A:1016592219341) algorithms with a sufficient [convergence and accuracy](https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be). Interested readers can also find the [stability](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://link.springer.com/article/10.1007/s00607-003-0033-3) fractional differential equations.

Let us suppose the following initial value problem with the Caputo fractional derivative <img src="https://latex.codecogs.com/svg.image?{}_{C}\!D_{t_0}^\beta" title="{}_{C}\!D_{t_0}^\beta" /> when <img src="https://latex.codecogs.com/svg.image?\beta>0" title="\beta>0" />

<img src="https://latex.codecogs.com/svg.image?{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))" title="{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))" />

with the initial condition <img src="https://latex.codecogs.com/svg.image?y(t_0)=y_0,y^{(1)}(t_0)=y^{(1)}_0,...,y^{(m-1)}(t_0)=y^{(m-1)}_0" title="y(t_0)=y0" />, where m the upper integer of the order of derivative.

We solve the problem by using predector corrector method (the equation (14) from this [paper](https://www.mdpi.com/2227-7390/6/2/16#)).


## Installation
If Julia is installed correctly, you can import FdeSolver.jl as:

```julia
import Pkg; Pkg.add("FdeSolver")
```

## API

Example1:
[Fractional nonlinear equation]( https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be)
<img src="https://latex.codecogs.com/gif.latex?\footnotesize{{}_{C}\!D_{t_0}^{\beta}y(t)=\frac{40320}{\Gamma(9-\beta)}t^{8-\beta}-3\frac{\Gamma(5+\beta/2)}{\Gamma(5-\beta/2)}t^{4-\beta/2}+\frac{9}{4}\Gamma(\beta+1)+\left(\frac{3}{2}t^{\beta/2}-t^4\right)^3-y(t)^{3/2}}" />, 
<img src="https://latex.codecogs.com/gif.latex?0<\beta\leq1" /> ,
subject to the initial condition <img src="https://latex.codecogs.com/gif.latex?y(0)=0" />.
The exact solution is
<img src="https://latex.codecogs.com/gif.latex?y(t)=t^8-3t^{4+\beta/2}+9/4t^\beta" />.

```julia
using FdeSolver
using Plots

## inputs
tSpan = [0, 1] # [intial time, final time]
y0 = 0 # intial value
β = 0.9 #order of the derivative

# Equation
F(t, n, β, y)=(40320 ./ gamma(9 - β) .* t[n] .^ (8 - β) .- 3 .* gamma(5 + β / 2)
           ./ gamma(5 - β / 2) .* t[n] .^ (4 - β / 2) .+ 9/4 * gamma(β + 1) .+
           (3/2 .* t[n] .^ (β / 2) .- t[n] .^ 4) .^ 3 .- y[n] .^ (3/2))

## Numerical solution
t, Yapp = FDEsolver(F, tSpan, y0, β)

#plot
plot(t, Yapp, linewidth = 5, title = "Solution of a 1D fractional IVP",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")
plot!(t, t -> (t.^8 - 3 * t .^ (4 + β / 2) + 9/4 * t.^β),
      lw = 3, ls = :dash, label = "Exact solution")
```


Example2: 
[Lotka-volterra-predator-prey](https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html)

```julia
using FdeSolver
using Plots

## inputs
tSpan = [0, 25] # [intial time, final time]
y0 = [34, 6] # initial values
β = [0.98, 0.99] # order of derivatives

## ODE model

par = [0.55, 0.028, 0.84, 0.026] # model parameters

function F(t, n, β, y, par)

    α1=par[1] #growth rate of the prey population
    β1=par[2] #rate of shrinkage relative to the product of the population sizes
    γ=par[3] #shrinkage rate of the predator population
    δ=par[4] #growth rate of the predator population as a factor of the product
             #of the population sizes

    u= y[n, 1] #population size of the prey species at time t[n]
    v= y[n, 2] #population size of the predator species at time t[n]

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
```

Example3:
[SIR model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology)

One application of using fractional calculus is taking into account effects of [memory](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.022409) in modeling including epidemic evolution.
```julia
using FdeSolver
using Plots

## inputs
I0 = 0.001 #intial value of infected
tSpan = [0, 100] # [intial time, final time]
y0 = [1 - I0, I0, 0] # initial values [S0,I0,R0]
α = [1, 1, 1] # order of derivatives
h = 0.1 # step size of computation (default=0.01)

## ODE model
par = [0.4, 0.04] # parameters [β, recovery rate]

function F(t, n, α, y, par)

    #parameters
    β = par[1] #infection rate
    γ = par[2] #recovery rate

    S = y[n, 1] #Susceptible
    I = y[n, 2] #Infectious
    R = y[n, 3] #Recovered

    #System equation
    dSdt = - β .* S .* I
    dIdt = β .* S .* I .- γ .* I
    dRdt = γ .* I

    return [dSdt, dIdt, dRdt]

end
## Solution and plotting
t, Yapp = FDEsolver(F, tSpan, y0, α, par, h = h)

plot(t, Yapp, linewidth = 5, title = "Numerical solution of SIR model",
     xaxis = "Time (t)", yaxis = "SIR populations", label=["Susceptible" "Infectious" "Recovered"])
```
