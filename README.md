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
[Fractional Riccati equation](https://www.sciencedirect.com/science/article/pii/S0898122111003245#s000055)
<img src="https://latex.codecogs.com/gif.latex?{}_{C}\!D_{t_0}^{\beta}y(t)=1+2y(t)-[y(t)]^2" />, <img src="https://latex.codecogs.com/gif.latex?0<\beta\leq1" /> ,
subject to the initial condition <img src="https://latex.codecogs.com/gif.latex?y(0)=0" />.
The exact solution for <img src="https://latex.codecogs.com/gif.latex?\beta=1" /> is

<img src="https://latex.codecogs.com/gif.latex?y(t)=1+\sqrt{2}tanh\left(\sqrt{2}t+\frac{1}{2}ln\Bigg(\frac{\sqrt{2}-1}{\sqrt{2}+1}\Bigg)\right)." />

```julia
using FdeSolver
using SpecialFunctions
using Plots

tSpan = [0, 1.5]
y0 = 0
β = 0.9

function F(t, n, β, y)

    return (40320 ./ gamma(9 - β) .* t[n] .^ (8 - β) .- 3 .* gamma(5 + β / 2) ./ gamma(5 - β / 2) .* t[n] .^ (4 - β / 2) .+ 9/4 * gamma(β + 1) .+ (3/2 .* t[n] .^ (β / 2) .- t[n] .^ 4) .^ 3 .- y[n] .^ (3/2))

end

t, Yapp = FDEsolver(F, tSpan, y0, β)

plot(t, Yapp, linewidth = 5, title = "Solution of a system of 2 FDEs", xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation 1")
plot!(t, t -> (t.^8 - 3 * t .^ (4 + β / 2) + 9/4 * t.^β), lw = 3, ls = :dash, label = "Exact solution 1")
```


Example2: 
```julia
using FdeSolver
using Plots
using SpecialFunctions

tSpan = [0, 5]
y0 = [1, 0.5, 0.3]
β = [0.5, 0.2, 0.6]

function F(t, n, β, y)

    F1 = 1 / sqrt(pi) * (((y[n, 2] - 0.5) * (y[n, 3] - 0.3))^(1 / 6) + t[n]^(1 / 2))
    F2 = gamma(2.2) * (y[n, 1] - 1)
    F3 = gamma(2.8) / gamma(2.2) * (y[n, 2] - 0.5)

    return [F1, F2, F3]

end

t, Yapp = FDEsolver(F, tSpan, y0, β)

plot(t, Yapp, linewidth = 5, title = "Solution of system 33",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")

plot!(t, t -> (t .+ 1), lw = 3, ls = :dash, label = "Exact solution")
plot!(t, t -> (t.^1.2 .+ 0.5), lw = 3, ls = :dash, label = "Exact solution")
plot!(t, t -> (t.^1.8 .+ 0.3), lw = 3, ls = :dash, label = "Exact solution")
```


Example3:
[Lotka-volterra-predator-prey](https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html)
```julia
using FdeSolver
using Plots

## inputs
tSpan = [0, 25] # [intial time, final time]
y0 = [34, 6] # initial values
β = [0.98, 0.99] # order of derivatives

## ODE model
α1=0.55 #growth rate of the prey population
β=0.028 #rate of shrinkage relative to the product of the population sizes
γ=0.84 #shrinkage rate of the predator population
δ=0.026 #growth rate of the predator population as a factor of the product
        #of the population sizes
        
par = [α1, β1, γ, δ] # model parameters

function F(t, n, β, y, par)

    F1 = par[1] .* y[n, 1] .- par[2] .* y[n, 1] .* y[n, 2]
    F2 = - par[3] .* y[n, 2] .+ par[4] .* y[n, 1] .* y[n, 2]

    [F1, F2]

end
## Solution of the 
t, Yapp = FDEsolver(F, tSpan, y0, β, par)

plot(t, Yapp, linewidth = 5, title = "Solution to LV model with 2 FDEs",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")

```

