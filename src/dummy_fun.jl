# this dummy function is meant to keep the second positional argument of _FDEsolver,
# the Jacobian matrix J, jidden from the user, so that it is not necessary to give
# nothing as the second argument when the Jacobian matrix is not provided.

function FDEsolver(F, tSpan, y0, β, par...; h = 2^-5, nc = 1, J = nothing, StopIt = "Standard", tol = 10e-6, itmax = 100)

    if isnothing(J)

        _FDEsolver(F, nothing, tSpan, y0, β, par..., h = h, nc = nc, StopIt = StopIt, tol = tol, itmax = itmax)

    else

        _FDEsolver(F, J, tSpan, y0, β, par..., h = h, nc = nc, StopIt = StopIt, tol = tol, itmax = itmax)

    end

end
