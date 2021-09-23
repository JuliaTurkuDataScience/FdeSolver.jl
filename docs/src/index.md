```@meta
CurrentModule = FdeSolver
```

# FdeSolver

Documentation for [FdeSolver](https://github.com/JuliaTurkuDataScience/FdeSolver.jl).

```@docs
greet()
FDEsolver(F, tSpan, y0, β, par...; h = 0.01, nc = 2, tol = 10^(-6), itmax = 30)
```

This is a Pkg in **Julia** for solution to a class of fractional differential equations and system equations.
Many advanced source codes are available in [MATLAB](https://www.dm.uniba.it/members/garrappa/software), but they are not open source projects like this one in Julia. Hence, the purpose is to develop a Julia package that numerically solves nonlinear fractional ordinary differential equations.

```@index
```

```@autodocs
Modules = [FdeSolver]
```
