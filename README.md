# FdeSolver.jl: Solving fractional differential equations

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaturkudatascience.github.io/FdeSolver.jl/stable/readme/)
[![CI](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl/branch/main/graph/badge.svg?token=SJ5F6RQ31P)](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7461985.svg)](https://doi.org/10.5281/zenodo.7461985)

This is a **Julia** package for fractional differential equations and ODEs. It provides numerical solutions for nonlinear fractional ordinary differential equations (in the sense of Caputo).

Related work includes the independent [FractionalDiffEq](https://github.com/SciFracX/FractionalDiffEq.jl) Julia Package that provides solutions of differential equations with different fractional operators, and earlier tools that are available in [Matlab](https://www.dm.uniba.it/members/garrappa/software).

## Publications

The package development is linked with the following publications/preprints:

- Three-species Lotka-Volterra model with respect to Caputo and Caputo-Fabrizio fractional operators. M Khalighi, L Eftekhari, S Hosseinpour, L Lahti. Symmetry 13 (3):368, 2021. [doi:10.3390/sym13030368](https://doi.org/10.3390/sym13030368)

- Quantifying the impact of ecological memory on the dynamics of interacting communities. M Khalighi, G Sommeria-Klein, D Gonze, K Faust, L Lahti. bioRxiv 2021.09.01.458486 [doi:10.1101/2021.09.01.458486](https://doi.org/10.1101/2021.09.01.458486)


## Method

We implement the [predictor-corrector](https://link.springer.com/article/10.1023/A:1016592219341) algorithms with a sufficient [convergence and accuracy](https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be), including fast Fourier transform technique that gives us high computation speed. Interested readers can also find the [stability](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://link.springer.com/article/10.1007/s00607-003-0033-3) fractional differential equations.

Let us suppose the following initial value problem with the Caputo fractional derivative ``{}_{C}\!D_{t_0}^\beta``  for  ``\beta>0``:

```math
{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))
```

with the initial condition:
```math
y(t_0)=y_0,y^{(1)}(t_0)=y^{(1)}_0,...,y^{(m-1)}(t_0)=y^{(m-1)}_0
```
where m is the smallest integer or equal to the order of derivative.

We solve the problem by using predictor corrector method (the equation (14) from this [paper](https://www.mdpi.com/2227-7390/6/2/16#)).


## Installation
If Julia is installed correctly, you can import FdeSolver.jl as:

```julia
import Pkg; Pkg.add("FdeSolver")
```

A few methods on its usage are explained in [Examples](https://juliaturkudatascience.github.io/FdeSolver.jl/stable/examples/).
