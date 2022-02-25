# API

*Let's solve some differential equations!*

## Package features
- Solve fractional calculus problems

## Function Documentation
```@docs
FDEsolver(F::Function, tSpan::Vector{<:Real}, y0::Union{Real, Vector{<:Real}, Matrix{<:Real}}, β::Union{Real, Vector{<:Real}}, par...; h = 2^-6, nc = 2, JF = nothing, StopIt = "Standard", tol = 10e-6, itmax = 100)
```

To access the manual of `FDEsolver` from the Julia REPL, type:
```julia
?FDESolver
```
