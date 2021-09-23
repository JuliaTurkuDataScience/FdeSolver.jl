"""
    FDEsolver(F, tSpan, y0, β, par...; h = 0.01, nc = 3, tol = 10^(-9), itmax = 10)

Solves fractional differential equations with a predictor-corrector approach.

    FDEsolver(F, tSpan, y0, β, J = JacobF, par...; h = 0.01, nc = 3, tol = 10^(-9), itmax = 10)

Additionally takes the Jacobian of the system (J) to evaluate the solution.

# Arguments
- `F`: the right side of the system of differential equations. It must be expressed
   in the form of a function and return a vector function with the same number of
   entries of order of derivatives. This function can also include a vector of
   parameters: par... .
- `tSpan::Vector{Number}`: the time span along which computation is performed.
- `y0`: the initial values in the form of a row vector (Vector{Number}) for β <= 1
   and a matrix (Matrix{Number}) for β > 1, where each column corresponds to the
   initial values of one differential equation and each row to the order of derivation.
- `β::Vector{Number}`: the orders of derivation in the form of a row vector, where
   each element corresponds to the order of one differential equation. It can take
   decimal as well as integer values.
- `J::Matrix{Number}`: the Jacobian of F.
- `par...`: additional parameters for the function F.
- `h::Number`: the step size for correction.
- `nc:Int64`: the desired number of corrections.
- `tol::Float64`: the tolerance.
- `ìtmax::Int64`: the maximal number of iterations.

# Examples
Solution of Predator-Prey Model with additional parameters as par:
```jldoctest
## Inputs
tSpan = [0, 25]     # [intial time, final time]
y0 = [34, 6]        # initial values
β = [0.98, 0.99]    # order of derivatives

## ODE Model

par = [0.55, 0.028, 0.84, 0.026] # model parameters

function F(t, n, β, y, par)

    α1 = par[1]      # growth rate of the prey population
    β1 = par[2]      # rate of shrinkage relative to the product of the population sizes
    γ = par[3]       # shrinkage rate of the predator population
    δ = par[4]       # growth rate of the predator population as a factor of the product
                     # of the population sizes

    u = y[n, 1]      # population size of the prey species at time t[n]
    v = y[n, 2]      # population size of the predator species at time t[n]

    F1 = α1 .* u .- β1 .* u .* v
    F2 = - γ .* v .+ δ .* u .* v

    [F1, F2]

end

## Solution
t, Yapp = FDEsolver(F, tSpan, y0, β, par)
```
Solution of SIR Model with and without Jacobian of the system:
```jldoctest
## Inputs
I0 = 0.001             # intial value of infected
tSpan = [0, 100]       # [intial time, final time]
y0 = [1 - I0, I0, 0]   # initial values [S0,I0,R0]
α = [1, 1, 1]          # order of derivatives
h = 0.1                # step size of computation (default=0.01)

## ODE Model
par = [0.4, 0.04]      # parameters [β, recovery rate]

function F(t, n, α, y, par)

    # Parameters
    β = par[1]         # infection rate
    γ = par[2]         # recovery rate

    S = y[n, 1]        # Susceptible
    I = y[n, 2]        # Infectious
    R = y[n, 3]        # Recovered

    # System Equations
    dSdt = - β .* S .* I
    dIdt = β .* S .* I .- γ .* I
    dRdt = γ .* I

    return [dSdt, dIdt, dRdt]

end

## Jacobian of ODE System
function JacobF(t, n, α, y, par)

    # Parameters
    β = par[1]         # infection rate
    γ = par[2]         # recovery rate

    S = y[n, 1]        # Susceptible
    I = y[n, 2]        # Infectious
    R = y[n, 3]        # Recovered

    # System Equations
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
t, Yapp = FDEsolver(F, tSpan, y0, α, par, h = h)
t1, Yapp1 = FDEsolver(F, tSpan, y0, α, J = JacobF, par, h = h)
```
"""
module FdeSolver

using SpecialFunctions
using LinearAlgebra

"""
    greet()
Motivates to get started with FdeSolver
"""
greet() = print("Hey, let's solve some FDEs!")

include("main.jl")
include("main_Jacob.jl")
include("SupFuns.jl")

export(FDEsolver)

end
