# FdeSolver

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaturkudatascience.github.io/FdeSolver.jl/stable/)
[![CI](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl/branch/main/graph/badge.svg?token=SJ5F6RQ31P)](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl)

This is a **Julia** package for fractional differential equations and ODEs. It provides numerical solutions for nonlinear fractional ordinary differential equations (in the sense of Caputo).

Related work includes the independent [FractionalCalculus](https://github.com/ErikQQY/FractionalCalculus.jl.git) Julia Package that covers fractional operators, and earlier tools that are available in [Matlab](https://www.dm.uniba.it/members/garrappa/software).

### Publications

The package development is linked with the following publications/preprints:

- Three-species Lotka-Volterra model with respect to Caputo and Caputo-Fabrizio fractional operators. M Khalighi, L Eftekhari, S Hosseinpour, L Lahti. Symmetry 13 (3):368, 2021. [doi:10.3390/sym13030368](https://doi.org/10.3390/sym13030368)

- Quantifying the impact of ecological memory on the dynamics of interacting communities. M Khalighi, D Gonze, K Faust, G Sommeria-Klein, L Lahti. bioRxiv 2021.09.01.458486 [doi:10.1101/2021.09.01.458486](https://doi.org/10.1101/2021.09.01.458486)


### Method

We implement the [predictor-corrector](https://link.springer.com/article/10.1023/A:1016592219341) algorithms with a sufficient [convergence and accuracy](https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be), including fast Fourier transform technique that gives us high computation speed. Interested readers can also find the [stability](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://link.springer.com/article/10.1007/s00607-003-0033-3) fractional differential equations.

Let us suppose the following initial value problem with the Caputo fractional derivative ${}_{C}\!D_{t_0}^\beta$ when $\beta>0$

<img src="https://latex.codecogs.com/svg.image?{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))" title="{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))" />

with the initial condition <img src="https://latex.codecogs.com/svg.image?y(t_0)=y_0,y^{(1)}(t_0)=y^{(1)}_0,...,y^{(m-1)}(t_0)=y^{(m-1)}_0" title="y(t_0)=y0" />, where m the smallest integer or equal to the order of derivative.

We solve the problem by using predictor corrector method (the equation (14) from this [paper](https://www.mdpi.com/2227-7390/6/2/16#)).


## Installation
If Julia is installed correctly, you can import FdeSolver.jl as:

```julia
import Pkg; Pkg.add("FdeSolver")
```

## API

Example1:
[Fractional nonlinear equation]( https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be)
<img src="https://latex.codecogs.com/gif.latex?0<\beta\leq1" /> ,
subject to the initial condition <img src="https://latex.codecogs.com/gif.latex?y(0)=0" />.
The exact solution is
<img src="https://latex.codecogs.com/gif.latex?y(t)=t^8-3t^{4+\beta/2}+9/4t^\beta" />.

```@repl 1
using FdeSolver
using Plots, SpecialFunctions

# Inputs
tSpan = [0, 1]     # [intial time, final time]
y0 = 0             # initial value
β = 0.9            # order of the derivative

# ODE Model
par = β
F(t, y, par) = (40320 ./ gamma(9 - par) .* t .^ (8 - par) .- 3 .* gamma(5 + par / 2)
           ./ gamma(5 - par / 2) .* t .^ (4 - par / 2) .+ 9/4 * gamma(par + 1) .+
           (3 / 2 .* t .^ (par / 2) .- t .^ 4) .^ 3 .- y .^ (3 / 2))

## Numerical solution
t, Yapp = FDEsolver(F, tSpan, y0, β, par)

# Plot
plot(t, Yapp, linewidth = 5, title = "Solution of a 1D fractional IVP",
     xaxis = "Time (t)", yaxis = "y(t)", label = "Approximation")
plot!(t, t -> (t.^8 - 3 * t .^ (4 + β / 2) + 9/4 * t.^β),
      lw = 3, ls = :dash, label = "Exact solution")
savefig("example1.png"); nothing # hide
```

![example1](example1.png)

Example2:
[Lotka-volterra-predator-prey](https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html)

```@repl 2
using FdeSolver
using Plots

# Inputs
tSpan = [0, 25]                    # [intial time, final time]
y0 = [34, 6]                       # initial values
β = [0.98, 0.99]                   # order of derivatives
par = [0.55, 0.028, 0.84, 0.026]   # model parameters

# ODE Model
function F(t, y, par)

    α1 = par[1]      # growth rate of the prey population
    β1 = par[2]      # rate of shrinkage relative to the product of the population sizes
    γ = par[3]       # shrinkage rate of the predator population
    δ = par[4]       # growth rate of the predator population as a factor of the product
                     # of the population sizes

    u = y[1]         # population size of the prey species at time t[n]
    v = y[2]         # population size of the predator species at time t[n]

    F1 = α1 .* u .- β1 .* u .* v
    F2 = - γ .* v .+ δ .* u .* v

    [F1, F2]

end

## Solution
t, Yapp = FDEsolver(F, tSpan, y0, β, par)

# Plot
plot(t, Yapp, linewidth = 5, title = "Solution to LV model with 2 FDEs",
     xaxis = "Time (t)", yaxis = "y(t)", label = ["Prey" "Predator"])
plot!(legendtitle = "Population of")
savefig("example2.png"); nothing # hide
```

![example2](example2.png)

Example3:
[SIR model](https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology)

One application of using fractional calculus is taking into account effects of [memory](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.022409) in modeling including epidemic evolution.

By defining the Jacobian matrix, the user can achieve a faster convergence based on the modified [Newton–Raphson](https://www.mdpi.com/2227-7390/6/2/16/htm) method.

```@repl 3
using FdeSolver
using Plots

# Inputs
I0 = 0.001             # intial value of infected
tSpan = [0, 100]       # [intial time, final time]
y0 = [1 - I0, I0, 0]   # initial values [S0,I0,R0]
α = [1, 1, 1]          # order of derivatives
h = 0.1                # step size of computation (default = 0.01)
par = [0.4, 0.04]      # parameters [β, recovery rate]

## ODE model
function F(t, y, par)

    # parameters
    β = par[1]    # infection rate
    γ = par[2]    # recovery rate

    S = y[1]   # Susceptible
    I = y[2]   # Infectious
    R = y[3]   # Recovered

    # System equation
    dSdt = - β .* S .* I
    dIdt = β .* S .* I .- γ .* I
    dRdt = γ .* I

    return [dSdt, dIdt, dRdt]

end

## Jacobian of ODE system
function JacobF(t, y, par)

    # parameters
    β = par[1]     # infection rate
    γ = par[2]     # recovery rate

    S = y[1]    # Susceptible
    I = y[2]    # Infectious
    R = y[3]    # Recovered

    # System equation
    J11 = - β * I
    J12 = - β * S
    J13 =  0
    J21 =  β * I
    J22 =  β * S - γ
    J23 =  0
    J31 =  0
    J32 =  γ
    J33 =  0

    J = [J11 J12 J13
         J21 J22 J23
         J31 J32 J33]

    return J

end

## Solution
t, Yapp = FDEsolver(F, tSpan, y0, α, par, JF = JacobF, h = h)

# Plot
plot(t, Yapp, linewidth = 5, title = "Numerical solution of SIR model",
     xaxis = "Time (t)", yaxis = "SIR populations", label=["Susceptible" "Infectious" "Recovered"])
savefig("example3.png"); nothing # hide
```

![example3](example3.png)

Example4:
Dynamics of interaction of N species microbial communities

The impact of [ecological memory](https://www.biorxiv.org/content/10.1101/2021.09.01.458486v1.abstract) on the dynamics of interacting communities can be quantified by solving fractional form ODE systems.

```@repl 4
using FdeSolver
using Plots

## inputs
tSpan = [0, 50]   # time span
h = 0.1           # time step
N = 20            # number of species
β = ones(N)       # order of derivatives
X0 = 2 * rand(N)  # initial abundances

# parametrisation
par = [2,
       2 * rand(N),
       rand(N),
       4 * rand(N, N),
       N]

# ODE model
function F(t, x, par)

    l = par[1] # Hill coefficient
    b = par[2] # growth rates
    k = par[3] # death rates
    K = par[4] # inhibition matrix
    N = par[5] # number of species

    Fun = zeros(N)

    for i in 1:N

        # inhibition functions
        f = prod(K[i, 1:end .!= i] .^ l ./
             (K[i, 1:end .!= i] .^ l .+ x[ 1:end .!= i] .^l))

        # System of equations
        Fun[i] = x[ i] .* (b[i] .* f .- k[i] .* x[ i])

    end

    return Fun

end

# Solution
t, Xapp = FDEsolver(F, tSpan, X0, β, par, h = h, nc = 3, tol = 10e-9)

# Plot
plot(t, Xapp, linewidth = 5,
     title = "Dynamics of microbial interaction model",
     xaxis = "Time (t)")
     yaxis!("Log Abundance", :log10, minorgrid = true)
savefig("example4.png"); nothing # hide
```

![example4](example4.png)
