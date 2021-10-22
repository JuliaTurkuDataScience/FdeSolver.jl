## some structures
struct JProblem
    ic
    f_fun
    problem_size::Int64
    param
    β # we should think about its type
    β_length::Int64
    JF
end

struct JMethod
    an::Matrix{Float64}
    a0::Matrix{Float64}
    hα1 # we should think about its type
    hα2 # we should think about its type
    μ::Int64
    μTol::Float64
    r::Int64
    StopIt::String
    itmax::Int64
end

struct JMethod_fft
    an_fft::Matrix{ComplexF64}
    index_fft::Matrix{Int64}
end

##
# Based on Roberto Garrappa's codes
function JDisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy, zn, N, METH, METH_fft, Probl)

    # I defined them bc they couldn't sum themselves in a loop (e.g. nxi += r)
    nxi = 0
    nxf = 0
    nyi = 0
    nyf = 0

    nxi = nx0
    nxf = nx0 + L * r - 1
    nyi = ny0
    nyf = ny0 + L * r - 1
    is = 1

    s_nxi = Int64.(zeros(L))
    s_nxf = Int64.(zeros(L))
    s_nyi = Int64.(zeros(L))
    s_nyf = Int64.(zeros(L))

    s_nxi[is] = nxi
    s_nxf[is] = nxf
    s_nyi[is] = nyi
    s_nyf[is] = nyf

    i_triangolo = 0
    stop = 0
    stop = false # We have to define stop (e.g. by zeros) before the loop to keep changes inside the loop!

    while !stop

        stop = ((nxi + r - 1 == nx0 + L * r - 1) || (nxi + r - 1 >= Nr - 1)) # It stops when current triangle ends at the end of the square

        zn = JQuadrato(nxi, nxf, nyi, nyf, fy, zn, N, METH, METH_fft, Probl)

        y, fy = JTriangolo(nxi, nxi + r - 1, t, y, fy, zn, N, METH, Probl)
        i_triangolo += 1

        if !stop

            if nxi + r - 1 == nxf   # The triangle ends where the square ends -> you level down

                i_Delta = Int64.(ff[i_triangolo])
                Delta = i_Delta * r

                nxi = s_nxf[is] + 1
                nxf = s_nxf[is]  + Delta
                nyi = s_nxf[is] - Delta + 1
                nyf = s_nxf[is]

                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf

            else # The triangle ends before the square -> a square is made next to it

                nxi += r
                nxf = nxi + r - 1
                nyi = nyf + 1
                nyf += r

                is += 1

                s_nxi[is] = nxi
                s_nxf[is] = nxf
                s_nyi[is] = nyi
                s_nyf[is] = nyf

            end

        end

    end

    return y, fy

end

##
# Based on Roberto Garrappa's codes
function JQuadrato(nxi, nxf, nyi, nyf, fy, zn, N, METH, METH_fft, Probl)

    coef_end = nxf - nyi + 1
    i_fft = Int64.(log2(coef_end / METH.r))
    funz_beg = nyi + 1
    funz_end = nyf + 1
    Nnxf = min(N, nxf)

    vett_funz = zeros(Probl.problem_size, coef_end)

    vett_funz[:, 1:coef_end] = ifelse(nyi == 0,
    [zeros(Probl.problem_size, 1) fy[:, funz_beg + 1:funz_end]], # Evaluation of the lowest square
    fy[:, funz_beg:funz_end])                                    # we have this "else" only in Jacob

    vett_funz_fft = fft(vett_funz, 2)
    zzn = zeros(Probl.problem_size, coef_end)

    for i in 1:Probl.problem_size

        i_β = min(Probl.β_length, i)
        Z = METH_fft.an_fft[i_β, METH_fft.index_fft[1, i_fft]:METH_fft.index_fft[2, i_fft]] .* vett_funz_fft[i, :]
        zzn[i, :] = real(ifft(Z))

    end

    zzn = zzn[:, nxf - nyf + 1:end]
    zn[:, nxi + 1:Nnxf + 1] += zzn[:, 1:Nnxf - nxi + 1]

    return zn

end

##
# Based on Roberto Garrappa's codes
function JTriangolo(nxi, nxf, t, y, fy, zn, N, METH, Probl)

    for n in nxi:min(N, nxf)

        # Evaluation of the predictor
        Φ = zeros(Probl.problem_size, 1)

        for j in nxi:n - 1

            Φ += METH.an[1:Probl.β_length,n - j + 1] .* fy[:, j + 1]

        end

        St = taylor_expansion(t[n + 1], Probl.ic)

        Φ_n = St + METH.hα2 .* (METH.a0[1:Probl.β_length, n + 1] .* fy[:, 1] .+ zn[:, n + 1] .+ Φ)
        yn0 = y[:, n]
        fn0 = f_value(Probl.f_fun(t[n + 1], yn0, Probl.param...), Probl.problem_size)
        Jfn0 = Probl.JF(t[n + 1], yn0, Probl.param...)
        Gn0 = yn0 - METH.hα2 .* fn0 - Φ_n

        stop = 0
        stop = false
        mu_it = 0

        yn1 = zeros(Probl.problem_size, 1)
        fn1 = zeros(Probl.problem_size, 1)

        while !stop

            JGn0 = I(Probl.problem_size) .- Diagonal(METH.hα2) * Jfn0
            yn1 = yn0 - inv(JGn0) * Gn0
            fn1 = f_value(Probl.f_fun(t[n + 1], yn1, Probl.param...), Probl.problem_size)
            Gn1 = yn1 - METH.hα2 .* fn1 - Φ_n

            mu_it += 1

            if METH.StopIt == "Convergence"

                stop = (norm(yn1 - yn0, Inf) < METH.μTol || norm(Gn1, Inf) <  METH.μTol)

                if (mu_it > METH.itmax && !stop)

                    stop = true

                end

            else

                stop = (mu_it == METH.μ)

            end

            yn0 = yn1
            Gn0 = Gn1

            if !stop

                Jfn0 = Probl.JF(t[n + 1], yn0, Probl.param...)

            end

        end

        y[:, n + 1] = yn1
        fy[:, n + 1] = fn1

    end

    return y, fy

end
