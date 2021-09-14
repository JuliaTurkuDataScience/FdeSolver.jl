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

```julia>
import Pkg; Pkg.add("FdeSolver")
```

## API

```julia>
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
