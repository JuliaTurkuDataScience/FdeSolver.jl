# FDE_Solver
This is the first version of a solver function for a class of fractional differential equations.
There are some related source codes in [MATLAB](https://www.mdpi.com/2227-7390/6/2/16) but not yet in Julia. Hence, the purpose is to develop a Julia package that numerically solves nonlinear fractional ordinary differential equations.
We implement the predictor-corrector algorithms.
You can find the details of the methods [here](https://link.springer.com/article/10.1023/A:1016592219341) that the authors have discussed the [convergence and accuracy](https://link.springer.com/article/10.1023/B:NUMA.0000027736.85078.be).
Interested readers can also find the [stability](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) of the methods and see how to implement the methods for solving [multi-term](https://link.springer.com/article/10.1007/s00607-003-0033-3) fractional differential equations.

## Method
Let us suppose the following initial value problem with the Caputo fractional derivative ${}_{C}\!D_{t_0}^\beta$ when $0<\beta<1$

${}_{C}\!D_{t_0}^{\beta}y(t)=f(t,y(t))$

with the initial condition $y(t_0)=y0$.

We use the equation (8) from this [paper](https://www.tandfonline.com/doi/full/10.1080/00207160802624331) to solve this problem:

$\Phi_n=a_{n,0}f(t_0,y_0)+\sum_{j=1}^{n-1}\alpha_{n-j}f(t_j,y_j)$

$y_n^P=y_0+h^{\beta}(\Phi_n-{\alpha}_{0}f(t_{n-2},y_{n-2})+2{\alpha_0}f(t_{n-1},y_{n-1}))$

$y_n=y_0+h^{\beta}(\Phi_n+{\alpha_0}f(t_{n},y_{n}^P))$

where $$\alpha_n=\frac{1}{\Gamma(\beta+2)}$$ if $n=0$, and $$\alpha_n=\frac{(n-1)^{\beta+1}-2n^{\beta+1}+(n+1)^{\beta+1}}{\Gamma(\beta+2)}$$ if $n= 1, 2, ...$, and

$a_{n,0}=\frac{(n-1)^{\beta+1}-n^{\beta}(n-\beta-1)}{\Gamma(\beta+2)}$


This is a 2-step method requiring two starting values, $y_0$ and $y_1$.
As we have the intial value $y_0$, only $y_1$ is needed.
Based on [equations (3) and (4)](https://www.tandfonline.com/doi/full/10.1080/00207160802624331), it can be calculated in 2 steps,
predictor:

$y_1^P=y_0+h^{β}\frac{f(t_0,y_0)}{\Gamma(\beta+1)}$

corrector:

$y_1=y_0+βh^{β}\frac{f(t_0,y_0)}{\Gamma(\beta+2)}+h^{β}\frac{f(t_0,y_1^P)}{\Gamma(\beta+2)}.$


## Usage
The solver is named FDESolver, and its arguments are listed below:

FDESolver(β=order of derivative,y0=initial value,t0=initial time,T=final time,h=time step,F=function $f(t,y)$ )

All arguments β, y0, t0, T, h should be a scaler, but the type of F should be a function.

output(t, y):  t is the discret-time, and y is the solution corresponding to the time.
