taylor_expansion(t0, t, y0) = y0

taylor_expansion(t0, t, y0::Matrix{Number}) = sum([(t - tSpan[1]) ^ k / factorial(k) * y0[:, k + 1] for k in 0:m - 1])
