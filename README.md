# FdeSolver.jl: Solving fractional differential equations

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaturkudatascience.github.io/FdeSolver.jl/stable/readme/)
[![CI](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl/branch/main/graph/badge.svg?token=SJ5F6RQ31P)](https://codecov.io/gh/JuliaTurkuDataScience/FdeSolver.jl)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7462094.svg)](https://doi.org/10.5281/zenodo.7462094)

This is a **Julia** package for fractional differential equations and ODEs. It provides numerical solutions for nonlinear fractional ordinary differential equations (in the sense of Caputo).

Related work includes the independent [FractionalDiffEq](https://github.com/SciFracX/FractionalDiffEq.jl) Julia Package that provides solutions of differential equations with different fractional operators, and earlier tools that are available in [Matlab](https://www.dm.uniba.it/members/garrappa/software).

The implemented models are generic and broadly applicable to modeling
multivariate signals from a single source or collected across multiple
sources. The dynamical models implemented in this package were
initially developed for modeling dynamics of interacting microbial
communities (Khalighi et al. 2022) but the models are more broadly
applicable and applicable to studying multi-omic and host-microbiome
interactions.



## Method

We implement the [predictor-corrector](https://doi.org/10.1023/A:1016592219341) algorithms with a sufficient [convergence and accuracy](https://doi.org/10.1023/B:NUMA.0000027736.85078.be), including fast Fourier transform technique that gives us high computation speed. Interested readers can also find the [stability](https://doi.org/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://doi.org/10.1007/s00607-003-0033-3) fractional differential equations.

Let us suppose the following initial value problem with the Caputo fractional derivative for  $\beta>0$:

```math
{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))
```

with the initial condition:
```math
y(t_0)=y_0,y^{(1)}(t_0)=y^{(1)}_0,...,y^{(m-1)}(t_0)=y^{(m-1)}_0
```
where m is the smallest integer or equal to the order of derivative.

We solve the problem by using [predictor corrector and Newton Raphson method](https://www.mdpi.com/2227-7390/6/2/16#).


## Installation
If Julia is installed correctly, you can import FdeSolver.jl as:

```julia
import Pkg; Pkg.add("FdeSolver")
```

A few methods on its usage are explained in [Examples](https://juliaturkudatascience.github.io/FdeSolver.jl/stable/examples/).


## Acknowledgments

We are grateful to all [contributors](https://github.com/JuliaTurkuDataScience/FdeSolver.jl/graphs/contributors). New issues and pull requests are welcome.

This research has received funding from
 * the Horizon 2020 Programme of the European Union within the [FindingPheno project](https://www.findingpheno.eu/) under grant agreement No 952914.
 * Research Council of Finland (grant 330887)

 Package documentation is compiled according to the guidelines provided in [PkgTutorial.jl](https://juliaturkudatascience.github.io/PkgTutorial.jl/dev/).

## Publications

**Kindly cite this work** as follows:

- Fdesolver: A Julia package for solving fractional differential equations M Khalighi, G Benedetti, L Lahti arXiv preprint, 2022.  [arXiv:2212.12550](https://arxiv.org/abs/2212.12550)


The package development is further linked with the following publications/preprints:

- Quantifying the impact of ecological memory on the dynamics of interacting communities. M Khalighi, G Sommeria-Klein, D Gonze, K Faust, L Lahti. PLoS Computational Biology 18(6), 2022 [doi:10.1371/journal.pcbi.1009396](https://doi.org/10.1371/journal.pcbi.1009396)

- Three-species Lotka-Volterra model with respect to Caputo and Caputo-Fabrizio fractional operators. M Khalighi, L Eftekhari, S Hosseinpour, L Lahti. Symmetry 13 (3):368, 2021. [doi:10.3390/sym13030368](https://doi.org/10.3390/sym13030368)

- Ebola epidemic model with dynamic population and memory, F Ndaïrou, M Khalighi, and L Lahti, Chaos, Solitons \& Fractals, 170: 113361, 2023. [doi:10.1016/j.chaos.2023.113361](https://doi.org/10.1016/j.chaos.2023.113361)






