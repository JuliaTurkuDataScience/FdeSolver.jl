using Revise
push!(LOAD_PATH, "./src")
using FdeSolver
using Plots

# time span
tSpan = [0, 100]

# time step
h = 0.1

# order of derivative
β = [1, 1, 1]

# growth rates
b = [1, 0.95, 1.05]

# death rates
k = [1, 1, 1]

# initial abundances
X = [0.1, 0.2, 0.3]

# parametrisation
par = [b[1] b[2] b[3];
       k[1] k[2] k[3]]

# system definition
function F(t, n, β, y, par)

    # Hill coefficient
    l = 2

    # inhibition matrix
    K = [nothing 0.1 0.1;
         0.1 nothing 0.1;
         0.1 0.1 nothing]

    # inhibition functions
    f1 = K[1, 2]^l / (K[1, 2]^l + y[n, 2]^l) * K[1, 3]^l / (K[1, 3]^l + y[n, 3]^l)
    f2 = K[2, 1]^l / (K[2, 1]^l + y[n, 1]^l) * K[2, 3]^l / (K[2, 3]^l + y[n, 3]^l)
    f3 = K[3, 1]^l / (K[3, 1]^l + y[n, 1]^l) * K[3, 2]^l / (K[3, 2]^l + y[n, 2]^l)

    # system of equations
    F1 = y[n, 1] * (par[1, 1] * f1 - par[2, 1] * y[n, 1])
    F2 = y[n, 2] * (par[1, 2] * f2 - par[2, 2] * y[n, 2])
    F3 = y[n, 3] * (par[1, 3] * f3 - par[2, 3] * y[n, 3])

    F = [F1, F2, F3]

end

# numerical solution
t, Yapp = FDEsolver(F, tSpan, X, β, par, h = h)

# plot
plot(t, Yapp, linewidth = 5,
     title = "3-member bacterial community with reciprocal inhibition",
     xaxis = "Time (t)", yaxis = "Relative Abundance (%)",
     label = ["X1" "X2" "X3"], legend = :topleft)
