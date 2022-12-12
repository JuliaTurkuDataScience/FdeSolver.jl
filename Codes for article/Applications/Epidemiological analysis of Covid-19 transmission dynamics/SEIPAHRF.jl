# ref:
# https://doi.org/10.1016/j.chaos.2020.109846
# https://doi.org/10.1016/j.chaos.2021.110652
# https://doi.org/10.3390/axioms10030135

using Optim, StatsBase
using FdeSolver
using Plots
using SpecialFunctions
using CSV
using DataFrames
using Dates

# Dataset
Death=CSV.read("time_series_covid19_deaths_global.csv", DataFrame)
Recover=CSV.read("time_series_covid19_recovered_global.csv", DataFrame)
Confirmed=CSV.read("time_series_covid19_confirmed_global.csv", DataFrame)

D=Death[Death.Country .== "Spain",39:120]
R=Recover[Recover.Country .== "Spain",39:120]
C=Confirmed[Confirmed.Country .== "Spain",39:120]

I=C .- D .-R
data=hcat(Vector(I[1,:]),Vector(R[1,:]),Vector(D[1,:]))


C=Float64.([6, 12, 19, 25, 31, 38, 44, 60, 80, 131, 131,259, 467, 688, 776,
 1776, 1460, 1739, 1984, 2101,2590, 2827, 3233, 3892, 3697, 3151,
 3387, 2653, 2984, 2473, 2022,1820, 1998, 1506, 1278, 2051, 1772,
 1891, 399, 894, 397, 650, 415,518, 412, 439, 441, 435, 579, 206,
 130, 120, 143, 146, 102, 46, 45,20, 31, 26, 11, 18, 27, 29, 39, 39])
D=[0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8, 15, 15, 25, 26, 26, 38, 43, 46,
 45, 57, 64, 66, 73, 73, 86, 89, 97, 108, 97, 254, 121, 121, 142, 106,
 106, 98, 115, 118, 109, 97, 150, 71, 52, 29, 44, 37, 35, 42, 31, 38,
 31, 30, 28, 27, 23, 17, 22, 11, 7, 14, 10, 14, 13, 13]
## System definition

# Define SIR model.
function SIR(t, u, p)
    # Model parameters.
	β=2.55
	l=1.56
	β′=7.65
	κ=0.25
	ρ₁=0.58
	ρ₂=0.001
	γₐ=0.94
	γᵢ=0.27
	γᵣ=0.5
	δᵢ=1/23
	δₚ=1/23
	δₕ=1/23

	p=N

    # Current state.
    S, E, I, P, A, H, R, F = u

# ODE
    dS = - β * I * S/N - l * β * H * S/N - β′* P * S/N
    dE = β * I * S/N + l * β * H * S/N + β′ *P* S/N - κ * E
    dI = κ * ρ₁ * E - (γₐ + γᵢ )*I - δᵢ * I
    dP = κ* ρ₂ * E - (γₐ + γᵢ)*P - δₚ * P
    dA = κ *(1 - ρ₁ - ρ₂ )* E
	dH = γₐ *(I + P ) - γᵣ *H - δₕ *H
	dR = γᵢ * (I + P ) + γᵣ* H
	dF = δᵢ * I + δₚ* P + δₕ *H
    return [dS, dE, dI, dP, dA, dH, dR, dF]
end

#initial conditions
N=11000000/250
par=N
S0=N-6; E0=0; I0=1; P0=5; A0=0; H0=0; R0=0; F0=0
X0=[S0, E0, I0, P0, A0, H0, R0, F0]
tspan=[1,66]
t1, x1 = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
plot(t1,x1)

XX=sum(x1[1:10:end,[3,4,6]], dims=2)
plot(t1[1:10:end],XX)
scatter!(C)

function loss(μ)
	if size(X0,2) != Int64(ceil(maximum(μ)))
		indx=findall(x-> x>1, μ)
		μ[indx]=ones(length(indx))
	end
	_, x = FDEsolver(SIR, tspan, X0, μ, par, h = .1)
    appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    rmsd(C, appX; normalize=:true)
    # sqrt(sum(abs2, appX .- data1))
end

p_lo=.5*ones(8)
p_up=ones(8)
pvec=.9*ones(8)
Result=optimize(loss,p_lo,p_up,pvec,Fminbox(LBFGS()),
# Result=optimize(loss,p_lo,p_up,pvec,SAMIN(rt=.99),
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
pp=vcat(Optim.minimizer(Result))
pp[3]=1;pp[7]=1
tf, xf = FDEsolver(SIR, tspan, X0, pp, par, h = .1)

Xf=sum(xf[1:10:end,[3,4,6]], dims=2)
plot(t1[1:10:end],Xf)
plot!(t1[1:10:end],XX)
scatter!(C)

Err1=rmsd(C, vec(XX); normalize=:true)
Errf=rmsd(C, vec(Xf); normalize=:true)
#
# julia> pp
# 8-element Vector{Float64}:
#  0.9590130843157745
#  0.9103147682079581
#  0.9999999977774652
#  0.5000000007931332
#  0.7499999564104551
#  0.999999995960407
#  0.7499999564104551
#  0.7499999564104551
