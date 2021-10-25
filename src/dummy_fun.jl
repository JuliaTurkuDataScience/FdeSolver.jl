# this dummy function is meant to keep the second positional argument of _FDEsolver,
# the Jacobian matrix J, jidden from the user, so that it is not necessary to give
# nothing as the second argument when the Jacobian matrix is not provided.

function FDEsolver(F, tSpan, y0, β, par...; h = 2^-6, nc = 2, J = nothing, StopIt = "Standard", tol = 10e-6, itmax = 100)

    pos_args = PositionalArguments(F, J, tSpan, y0, β)
    opt_args = OptionalArguments(h, nc, StopIt, tol, itmax)

    _FDEsolver(pos_args.F, pos_args.JF, pos_args.tSpan, pos_args.y0, pos_args.β, par..., h = opt_args.h, nc = opt_args.nc, StopIt = opt_args.StopIt, tol = opt_args.tol, itmax = opt_args.itmax)

end
