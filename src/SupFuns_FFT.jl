## define initial values ##
# I changed this part. We defineY0 only for intial conditions that only used in StartingTerm (taylor_expansion)
function defineY0(y0, β)

    Y0 = zeros(size(y0,1),Int64.(ceil(maximum(β))))

    Y0[:,1] .= y0

end

function defineY0(y0::Vector{<:Real}, β)

    Y0 = zeros(size(y0,1),Int64.(ceil(maximum(β))))

    Y0[:,1] .= y0
end

function defineY0(y0::Matrix{<:Real}, β) # I am not sure about it yet ks well if it works well for order>1

    Y0 = zeros(size(y0,1),Int64.(ceil(maximum(β))))
    Y0[: , 1] .= y0[:,1]

end

## Gamma function for vectors ##
# Γ(b) = map(gamma, b)

##
# Based on Roberto Garrappa's codes
function DisegnaBlocchi(L, ff, r, Nr, nx0, ny0, t, y, fy,zn_pred, zn_corr, N , METH, METH_fft, Probl)

nxi = zeros(1); nxf = zeros(1); nyi = zeros(1); nyf = zeros(1) # I defined them bc they couldn't sum themselves in a loop (e.g. nxi += r)
nxi = nx0 ; nxf = nx0 + L*r - 1
nyi = ny0 ; nyf = ny0 + L*r - 1
is = 1
s_nxi= Int64.(zeros(L)); s_nxf = Int64.(zeros(L)); s_nyi = Int64.(zeros(L)); s_nyf = Int64.(zeros(L))
s_nxi[is] = nxi ; s_nxf[is] = nxf ; s_nyi[is] = nyi ; s_nyf[is] = nyf

i_triangolo = 0 ;
stop = zeros(1); stop = false # We have to define stop (e.g. by zeros) before the loop to keep changes inside the loop!
while ~stop

    stop = (nxi+r-1 == nx0+L*r-1) || (nxi+r-1>=Nr-1) # It stops when current triangle ends at the end of the square

    zn_pred, zn_corr = Quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, N, METH, METH_fft, Probl)

    y, fy = Triangolo(nxi, nxi+r-1, t, y, fy, zn_pred, zn_corr, N, METH, METH_fft, Probl)
    i_triangolo += 1

    if ~stop
        if nxi+r-1 == nxf   # The triangle ends where the square ends -> you level down
            i_Delta = Int64.(ff[i_triangolo])
            Delta = i_Delta*r
            nxi = s_nxf[is]+1 ; nxf = s_nxf[is]  + Delta
            nyi = s_nxf[is] - Delta +1; nyf = s_nxf[is]
            s_nxi[is] = nxi ; s_nxf[is] = nxf ; s_nyi[is] = nyi ; s_nyf[is] = nyf
        else # The triangle ends before the square -> a square is made next to it
            nxi += r ; nxf = nxi + r - 1 ; nyi = nyf + 1 ; nyf += r
            is += 1
            s_nxi[is] = nxi ; s_nxf[is] = nxf ; s_nyi[is] = nyi ; s_nyf[is] = nyf
        end
    end

end

return y, fy
end

##
# Based on Roberto Garrappa's codes
function Quadrato(nxi, nxf, nyi, nyf, fy, zn_pred, zn_corr, N, METH, METH_fft, Probl)

coef_end = nxf-nyi+1
i_fft = Int64.(log2(coef_end/METH.r))
funz_beg = nyi+1 ; funz_end = nyf+1
Nnxf = min(N,nxf)

Indxfft=Int64.(METH_fft.index_fft) # IDK why METH_fft.index_fft was truned to complex number! So we need to change the type to Int.
 # Evaluation convolution segment for the predictor
vett_funz=zeros(Probl.problem_size,coef_end)
vett_funz[:,1:funz_end-funz_beg+1] = fy[:,funz_beg:funz_end]
vett_funz_fft = fft(vett_funz,2)
zzn_pred = zeros(Probl.problem_size,coef_end)
for i = 1 : Probl.problem_size
    i_β = min(Probl.β_length,i)
    Z = METH_fft.bn_fft[i_β,Indxfft[1,i_fft]:Indxfft[2,i_fft]].*vett_funz_fft[i,:]
    zzn_pred[i,:] = real(ifft(Z))
end
zzn_pred = zzn_pred[:,nxf-nyf:end-1]
zn_pred[:,nxi+1:Nnxf+1] += zzn_pred[:,1:Nnxf-nxi+1]

 # Evaluation convolution segment for the corrector
if METH.μ > 0
    if nyi == 0 # Evaluation of the lowest square
        vett_funz = zeros(Probl.problem_size,coef_end)
        vett_funz[:, 1:funz_end-funz_beg+1] = [zeros(Probl.problem_size,1) fy[:,funz_beg+1:funz_end] ]
        vett_funz_fft = fft(vett_funz,2)
    end
    zzn_corr = zeros(Probl.problem_size,coef_end)
    for i = 1 : Probl.problem_size
        i_β = min(Probl.β_length,i)
        Z = METH_fft.an_fft[i_β,Indxfft[1,i_fft]:Indxfft[2,i_fft]].*vett_funz_fft[i,:]
        zzn_corr[i,:] = real(ifft(Z))
    end
    zzn_corr = zzn_corr[:,nxf-nyf+1:end]
    zn_corr[:,nxi+1:Nnxf+1] += zzn_corr[:,1:Nnxf-nxi+1]
else
    zn_corr = 0
end

return zn_pred, zn_corr

end

##
# Based on Roberto Garrappa's codes
function Triangolo(nxi, nxf, t, y, fy, zn_pred, zn_corr, N, METH, METH_fft, Probl)

for n in nxi:min(N,nxf)

    # Evaluation of the predictor
    Φ = zeros(Probl.problem_size,1)
    if nxi == 1 # Case of the first triangle
        j_beg = 0
    else # Case of any triangle but not the first
        j_beg = nxi
    end
    for j in j_beg:n-1
        Φ += METH.bn[1:Probl.β_length,n-j].*fy[:,j+1]
    end
    St = StartingTerm(t[n+1], Probl.ic)
    y_pred = St .+ METH.hα1.*(zn_pred[:,n+1] .+ Φ)
    f_pred = Probl.f_fun(t, y_pred, Probl.param...)

    # Evaluation of the corrector
    if METH.μ == 0
        y[:,n+1] = y_pred
        fy[:,n+1] = f_pred
    else
        j_beg = nxi
        Φ = zeros(Probl.problem_size,1) ;
        for j in j_beg:n-1
            Φ += METH.an[1:Probl.β_length,n-j+1].*fy[:,j+1]
        end
        Φ_n = St +
            METH.hα2.*(METH.a0[1:Probl.β_length,n+1].*fy[:,1] .+ zn_corr[:,n+1] .+ Φ)
        yn0 = y_pred ; fn0 = f_pred
        stop = zeros(1); stop = false; # it is defined for counting in the following loop while
        mu_it = 0
        yn1 = zeros(Probl.problem_size,1) ; fn1 = zeros(Probl.problem_size,1)
        while ~stop
            yn1 = Φ_n + METH.hα2.*fn0
            mu_it += 1
            if METH.StopIt == "Convergence"
                stop = norm(yn1-yn0,Inf) < METH.μTol
                if mu_it > METH.itmax && ~stop
                    stop = true
                end
            else
                stop = mu_it == METH.μ
            end
            fn1 = Probl.f_fun(t, yn1, Probl.param...)
            yn0 = yn1 ; fn0 = fn1
        end
        y[:,n+1] = yn1
        fy[:,n+1] = fn1
      end
end

return y, fy
end

##
# Based on Roberto Garrappa's code (if we remove our previous code without FFT, then we can changed this name to taylor_expansion)
function StartingTerm(t, ic)

ys = zeros(size(ic.y0,1),1)
for k in 1:maximum(ic.m_β)
    if length(ic.m_β) == 1
        ys = ys .+ (t-ic.t0).^(k-1)./ic.m_β_factorial[k].*ic.y0[:,k]
    else
        i_β = findall(k .<= ic.m_β)
        ys[i_β,1] = ys[i_β,1] + (t-ic.t0)^(k-1)*ic.y0[i_β,k]./ic.m_β_factorial[i_β,k]
    end
end

return ys
end
