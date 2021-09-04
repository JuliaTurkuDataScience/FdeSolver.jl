defineY(N, y0) = zeros(N + 1, length(y0))

defineY(N, y0::Matrix{Number}) = zeros(N + 1, size(y0)[2])
