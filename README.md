# FdeSolver

This is a solver Pkg in Julia for a class of fractional differential equations.
There are some related source codes in [MATLAB](https://www.dm.uniba.it/members/garrappa/software) but not yet in Julia. Hence, the purpose is to develop a Julia package that numerically solves nonlinear fractional ordinary differential equations.
We implement the predictor-corrector algorithms.
You can find the details of the methods [here](https://link.springer.com/article/10.1023/A:1016592219341) that the authors have discussed the [convergence and accuracy](https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be).
Interested readers can also find the [stability](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://link.springer.com/article/10.1007/s00607-003-0033-3) fractional differential equations.

## Method
Let us suppose the following initial value problem with the Caputo fractional derivative <img src="https://latex.codecogs.com/svg.image?{}_{C}\!D_{t_0}^\beta" title="{}_{C}\!D_{t_0}^\beta" /> when <img src="https://latex.codecogs.com/svg.image?\beta>0" title="\beta>0" />

<img src="https://latex.codecogs.com/svg.image?{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))" title="{}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))" />

with the initial condition <img src="https://latex.codecogs.com/svg.image?y(t_0)=y_0,y^{(1)}(t_0)=y^{(1)}_0,...,y^{(m-1)}(t_0)=y^{(m-1)}_0" title="y(t_0)=y0" />, where m the upper integer of the order of derivative.

We solve the problem by using predector corrector method (the equation (14) from this [paper](https://www.mdpi.com/2227-7390/6/2/16#)).


## Usage
**We will later update this part**

The solver is named FDESolver, and its arguments are listed below:

FDESolver(β=order of derivative,y0=initial value,t0=initial time,T=final time,h=time step,F=function $f(t,y)$ )

All arguments β, y0, t0, T, h should be a scaler, but the type of F should be a function.

output(t, y):  t is the discret-time, and y is the solution corresponding to the time.
