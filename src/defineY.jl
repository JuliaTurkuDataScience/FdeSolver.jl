function defineY(N, y0)

    Y = zeros(N + 1, length(y0))
    Y[1, :] .= y0

    Y

end

function defineY(N, y0::Matrix{Number})

    Y = zeros(N + 1, size(y0)[2])
    Y[1, :] .= y0[1, :]

    Y

end
