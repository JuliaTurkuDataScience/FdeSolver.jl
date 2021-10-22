function _FDEsolver(F, ::Nothing, tSpan, y0, β, par...; h = 2^-6, nc = 1, StopIt = "Standard", tol = 10e-6, itmax = 100)

    # Check compatibility size of the problem with number of fractional orders
    # if size(y0, 1) != length(β)
    #     error("Size Is Not Compatible
    #                 Size %i of the problem as obtained from initial conditions
    #                 (i.e. the number of rows of Y0) not compatible with the
    #                 number %i of fractional orders for multi-order systems.
    #                 ", size(y0, 1), length(β))
    # end

    y0 = defineY0(y0, β)
    β_length = length(β)
    problem_size = size(y0, 1)

    if β_length ==1
        β = β * ones(problem_size, 1)
        β_length = problem_size

    end

    # Storage of initial conditions
    ic = initial_conditions(tSpan[1], y0, Int64.(map(ceil, β)), zeros(β_length, Int64.(ceil(maximum(β)))))

    for i in 1:β_length

        for j in 0:ic.m_β[i] - 1

            ic.m_β_factorial[i, j + 1] = factorial(j)

        end

    end

    # Storage of information on the problem
    Probl = Problem(ic, F, problem_size,par, β, β_length)

    # Time discretization
    N = Int64.(cld(tSpan[2] - tSpan[1], h))
    t = tSpan[1] .+ collect(0:N) .* h

    # Enter the initial values
    # Y = defineY(N, y0, β)

    # Check compatibility size of the problem with size of the vector field
    f_temp = f_value(F(t[1], y0[:,1], par...), Probl.problem_size)

    # Number of points in which to evaluate weights and solution
    r = 16
    Nr::Int64 = ceil((N + 1) / r) * r
    Qr::Int64 = ceil(log2(Nr / r)) - 1
    NNr::Int64 = 2^(Qr + 1) * r

    # Preallocation of some variables
    y = zeros(Probl.problem_size, N + 1)
    fy = zeros(Probl.problem_size, N + 1)
    zn_pred = zeros(Probl.problem_size, NNr + 1)

    zn_corr = ifelse(nc > 0, zeros(Probl.problem_size, NNr + 1), 0)

    # Evaluation of coefficients of the PECE method
    nvett = 0:NNr + 1
    bn = zeros(Probl.β_length, NNr + 1)
    an = zeros(Probl.β_length, NNr + 1)
    a0 = zeros(Probl.β_length, NNr + 1)

    for i_β in 1:Probl.β_length

        find_β = findall(β[i_β] == β[1:i_β - 1])

        if !isempty(find_β) # it is for speeding up the computations; we can use multilpe distpach

            bn[i_β, :] = bn[find_β[1], :]
            an[i_β, :] = an[find_β[1], :]
            a0[i_β, :] = a0[find_β[1], :]

        else

            nβ = nvett.^β[i_β]
            nβ1 = nβ .* nvett

            bn[i_β, :] = nβ[2:end] - nβ[1:end - 1]
            an[i_β, :] = [1; (nβ1[1:end - 2] - 2 * nβ1[2:end - 1] + nβ1[3:end])]
            a0[i_β, :] = [0; (nβ1[1:end-2] - nβ[2:end-1].*(nvett[2:end-1] .- β[i_β] .- 1))]

        end

    end

    METH = Method(bn, an, a0, h.^β ./ Γ(β .+ 1), h .^ β ./ Γ(β .+ 2), nc, tol, r, StopIt, itmax)

    # Evaluation of FFT of coefficients of the PECE method
    if Qr >= 0

        index_fft = Int64.(zeros(2, Qr + 1)) # I have tried index_fft::Int64 = zeros(2,Qr+1) and I got an error for converting Type!

        for l in 1:Qr + 1

            if l == 1

                index_fft[1,l] = 1
                index_fft[2, l] = r * 2

            else

                index_fft[1, l] = index_fft[2, l - 1] + 1
                index_fft[2, l] = index_fft[2, l - 1] + 2^l * r

            end

        end

        bn_fft = ComplexF64.(zeros(Probl.β_length, index_fft[2, Qr + 1]))
        an_fft = ComplexF64.(zeros(Probl.β_length, index_fft[2, Qr + 1]))

        for l in 1:Qr + 1

            coef_end = 2^l * r

            for i_β in 1:Probl.β_length

                find_β = findall(β[i_β] == β[1:i_β - 1])

                if !isempty(find_β)

                    bn_fft[i_β, index_fft[1, l]:index_fft[2, l]] = bn_fft[find_β[1], index_fft[1, l]:index_fft[2, l]]
                    an_fft[i_β, index_fft[1, l]:index_fft[2, l]] = an_fft[find_β[1], index_fft[1, l]:index_fft[2, l]]

                else

                    bn_fft[i_β, index_fft[1, l]:index_fft[2, l]] = fft(METH.bn[i_β, 1:coef_end])
                    an_fft[i_β, index_fft[1, l]:index_fft[2, l]] = fft(METH.an[i_β, 1:coef_end])

                end

            end

        end

        # Method_fft = @SLVector (:bn_fft,:an_fft,:index_fft)
        METH_fft = Method_fft(bn_fft, an_fft, Int64.(index_fft))
    end

    # Initializing solution and proces of computation
    y[:, 1] = y0[:, 1]
    fy[:, 1] = f_temp
    y, fy = Triangolo(1, r - 1, t, y, fy, zn_pred, zn_corr, N, METH, Probl)

    # Main process of computation by means of the FFT algorithm
    ff = zeros(1, 2^(Qr + 2))
    ff[1:2] = [0, 2]
    card_ff = 2

    nx0 = 0
    ny0 = 0

    for qr in 0:Qr

        L = 2^qr
        y, fy = DisegnaBlocchi(L, ff, r, Nr, nx0 + L * r, ny0, t, y, fy, zn_pred, zn_corr, N, METH, METH_fft, Probl)
        ff[1:2 * card_ff] = [ff[1:card_ff]; ff[1:card_ff]]
        card_ff = 2 * card_ff
        ff[card_ff] = 4 * L

    end

     # Evaluation solution in T when T is not in the mesh
    if tSpan[end] < t[N + 1]

        c = [tSpan[end] - t(N)] / h
        t[N + 1] = tSpan[end]
        y[:, N + 1] = (1 - c) * y[:, N] + c * y[:, N + 1]

    end

    t = t[1:N + 1]
    y = transpose(y[:, 1:N + 1])

    return t, y

end
