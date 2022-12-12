# reference for model
# https://doi.org/10.1016/j.chaos.2020.109846
# https://doi.org/10.1016/j.chaos.2021.110652
# https://doi.org/10.3390/axioms10030135
using Optim, StatsBase
using FdeSolver
using Plots
using SpecialFunctions
using CSV, HTTP, DataFrames, Dates

# Dataset
repo=HTTP.get("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv") # dataset of Covid from CSSE
dataset_CC = CSV.read(repo.body, DataFrame) # all data of confirmed
Confirmed=dataset_CC[dataset_CC[!,2].=="Portugal",45:121] #comulative confirmed data of Portugal from 3/2/20 to 5/17/20
C=diff(Float64.(Vector(Confirmed[1,:])))# Daily new confirmed cases

#preprocessing (map negative values to zero and remove outliers)
₋Ind=findall(C.<0)
C[₋Ind].=0.0
outlier=findall(C.>1500)
C[outlier]=(C[outlier.-1]+C[outlier.-1])/2

## System definition

# parameters
β=2.55 # Transmission coeﬃcient from infected individuals
l=1.56 # Relative transmissibility of hospitalized patients
β′=7.65 # Transmission coeﬃcient due to super-spreaders
κ=0.25 # Rate at which exposed become infectious
ρ₁=0.58 # Rate at which exposed people become infected I
ρ₂=0.001 # Rate at which exposed people become super-spreaders
γₐ=0.94 # Rate of being hospitalized
γᵢ=0.27 # Recovery rate without being hospitalized
γᵣ=0.5 # Recovery rate of hospitalized patients
δᵢ=1/23 # Disease induced death rate due to infected class
δₚ=1/23 # Disease induced death rate due to super-spreaders
δₕ=1/23 # Disease induced death rate due to hospitalized class
# Define SIR model
function SIR(t, u, par)
    # Model parameters.
	N, β, l, β′, κ,	ρ₁,	ρ₂,	γₐ,	γᵢ,	γᵣ,	δᵢ,	δₚ, δₕ=par

    # Current state.
    S, E, I, P, A, H, R, F = u

# ODE
    dS = - β * I * S/N - l * β * H * S/N - β′* P * S/N # susceptible individuals
    dE = β * I * S/N + l * β * H * S/N + β′ *P* S/N - κ * E # exposed individuals
    dI = κ * ρ₁ * E - (γₐ + γᵢ )*I - δᵢ * I #symptomatic and infectious individuals
    dP = κ* ρ₂ * E - (γₐ + γᵢ)*P - δₚ * P # super-spreaders individuals
    dA = κ *(1 - ρ₁ - ρ₂ )* E # infectious but asymptomatic individuals
	dH = γₐ *(I + P ) - γᵣ *H - δₕ *H # hospitalized individuals
	dR = γᵢ * (I + P ) + γᵣ* H # recovery individuals
	dF = δᵢ * I + δₚ* P + δₕ *H # dead individuals
    return [dS, dE, dI, dP, dA, dH, dR, dF]
end

#initial conditions
N=10280000/875 # Population Size
S0=N-5; E0=0; I0=4; P0=1; A0=0; H0=0; R0=0; F0=0
X0=[S0, E0, I0, P0, A0, H0, R0, F0] # initial values
tspan=[1,length(C)] # time span [initial time, final time]

par=[N, β,	l,	β′,	κ,	ρ₁,	ρ₂,	γₐ,	γᵢ,	γᵣ,	δᵢ,	δₚ, δₕ] # parameters

## optimazation of β for integer order model

function loss_1(b)# loss function
	par[2]=b[1]
	_, x = FDEsolver(SIR, tspan, X0, ones(8), par, h = .1)
    appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    rmsd(C, appX; normalize=:true) # Normalized root-mean-square error
end

p_lo_1=[1.4] #lower bound for β
p_up_1=[4.0] # upper bound for β
p_vec_1=[2.5] #  initial guess for β
Res1=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,Fminbox(BFGS()),# Broyden–Fletcher–Goldfarb–Shanno algorithm
# Result=optimize(loss_1,p_lo_1,p_up_1,p_vec_1,SAMIN(rt=.99), # Simulated Annealing algorithm (sometimes it has better perfomance than (L-)BFGS)
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
p1=vcat(Optim.minimizer(Res1))
par1=par; par1[2]=p1[1]

## optimazation of β and order of commensurate fractional order model
function loss_F_1(pμ)
	par[2] = pμ[1] # infectivity rate
	μ = pμ[2] # order of derivatives

	_, x = FDEsolver(SIR, tspan, X0, μ*ones(8), par, h = .1)
    appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    rmsd(C, appX; normalize=:true)
end

p_lo_f_1=vcat(1.4,.5)
p_up_f_1=vcat(4,1)
p_vec_f_1=vcat(2.5,.9)
ResF1=optimize(loss_F_1,p_lo_f_1,p_up_f_1,p_vec_f_1,Fminbox(LBFGS()), # LBFGS is suitable for large scale problems
# Result=optimize(loss,p_lo,p_up,pvec,SAMIN(rt=.99),
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
pμ=vcat(Optim.minimizer(ResF1))
parf1=par; parf1[2]=pμ[1]; μ1=pμ[2]

## optimazation of β and order of incommensurate fractional order model
function loss_F_8(pμ)
	par[2] = pμ[1] # infectivity rate
	μ = pμ[2:9] # order of derivatives
	if size(X0,2) != Int64(ceil(maximum(μ))) # to prevent any errors regarding orders higher than 1
		indx=findall(x-> x>1, μ)
		μ[indx]=ones(length(indx))
	end
	_, x = FDEsolver(SIR, tspan, X0, μ, par, h = .1)
    appX=vec(sum(x[1:10:end,[3,4,6]], dims=2))
    rmsd(C, appX; normalize=:true)
end

p_lo=vcat(1.4,.5*ones(8))
p_up=vcat(4,ones(8))
pvec=vcat(2.5,.9*ones(8))
ResF8=optimize(loss_F_8,p_lo,p_up,pvec,Fminbox(LBFGS()), # LBFGS is suitable for large scale problems
# Result=optimize(loss,p_lo,p_up,pvec,SAMIN(rt=.99),
			Optim.Options(outer_iterations = 10,
						  iterations=10000,
						  show_trace=true,
						  show_every=1))
pp=vcat(Optim.minimizer(ResF8))
parf8=par; parf8[2]=pp[1]; μ8=pp[2:9]
#
# the optimized β and orders
# julia> pp
# 9-element Vector{Float64}:
#  2.841678343258684
#  0.764514142673873
#  0.7248610635678789
#  0.8758563369131535
#  0.9999999997148052
#  0.74983538386169
#  0.611714850156578
#  0.74983538386169
#  0.74983538386169

## plotting
DateTick=Date(2020,3,3):Day(1):Date(2020,5,17)
DateTick2= Dates.format.(DateTick, "d u")

t1, x1 = FDEsolver(SIR, tspan, X0, ones(8), par1, h = .1) # solve ode model
_, xf1 = FDEsolver(SIR, tspan, X0, μ1*ones(8), parf1, h = .1) # solve commensurate fode model
_, xf8 = FDEsolver(SIR, tspan, X0, μ8, parf8, h = .1) # solve incommensurate fode model

X1=sum(x1[1:10:end,[3,4,6]], dims=2)
Xf1=sum(xf1[1:10:end,[3,4,6]], dims=2)
Xf8=sum(xf8[1:10:end,[3,4,6]], dims=2)

# Err1=rmsd(C, vec(X1); normalize=:true) # NRMSE for ode model
# Errf1=rmsd(C, vec(Xf1); normalize=:true) # NRMSE for commensurate fode model
# Errf8=rmsd(C, vec(Xf8); normalize=:true) # NRMSE for incommensurate fode model

Err1=rmsd(C, vec(X1)) # RMSD for ode model
Errf1=rmsd(C, vec(Xf1)) # RMSD for commensurate fode model
Errf8=rmsd(C, vec(Xf8)) # RMSD for incommensurate fode model

plot(DateTick2,X1, ylabel="Daily new confirmed cases",
     label="M1",xrotation=rad2deg(pi/3), linestyle=:dashdot)
    plot!(Xf1, label="Mf1")
    plot!(Xf8,  label="Mf8", linestyle=:dash)
    scatter!(C, label= "Real data",legendposition=(.85,1),legend=:false,
	title = "(b) Portugal" , titleloc = :left, titlefont = font(10))
	plPortugal=bar!(["M1" "Mf1" "Mf8"],[Err1 Errf1 Errf8], ylabel="RMSD",
		legend=:false, bar_width=2,yguidefontsize=8,xtickfontsize=7,
    inset = (bbox(0.04, 0.08, 70px, 50px, :right)),
    subplot = 2,
    bg_inside = nothing)
