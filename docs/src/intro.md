# FdeSolver

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaturkudatascience.github.io/FdeSolver.jl/stable/intro/)
[![CI](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl/branch/main/graph/badge.svg?token=SJ5F6RQ31P)](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl)

This is a Pkg in **Julia** for solution to a class of fractional differential equations and system equations.
Many advanced source codes are available in [MATLAB](https://www.dm.uniba.it/members/garrappa/software), but they are not open source projects like this one in Julia. Hence, the purpose is to develop a Julia package that numerically solves nonlinear fractional ordinary differential equations.

## Method

We implement the [predictor-corrector](https://link.springer.com/article/10.1023/A:1016592219341) algorithms that the convergence and accuracy has been studied [here](https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be). Interested readers can also find the [stability](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://link.springer.com/article/10.1007/s00607-003-0033-3) fractional differential equations.

Let us suppose the following initial value problem

```math
{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t)),
```
with the initial condition
```math
y(t_0)=y_0,y^{(1)}(t_0)=y^{(1)}_0,...,y^{(m-1)}(t_0)=y^{(m-1)}_0,
```
in which, ```math
{}_{C}\!D_{t_0}^\beta,  \beta>0,
```
the Caputo fractional derivative and m the upper integer of the order of derivative.

We solve the problem by using predictor-corrector method (the equation (14) from this [paper](https://www.mdpi.com/2227-7390/6/2/16#)).


## Installation
If Julia is installed correctly, you can import FdeSolver.jl as:

```julia
import Pkg; Pkg.add("FdeSolver")
# or
]add FdeSolver
```
