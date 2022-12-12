using Plots, CSV, DataFrames
using StatsPlots
using EasyFit, GLM # for plotting linear line in the background
using LaTeXStrings # for math in plots
# theme(:ggplot2)
theme(:default)
# Data from Matlab
M_Ex1 = CSV.read("data_matlab/RndEx1.csv", DataFrame, header = 1)
M_Ex2 = CSV.read("data_matlab/RndEx2.csv", DataFrame, header = 1)
M_Ex3 = CSV.read("data_matlab/RndEx3.csv", DataFrame, header = 1)
M_Ex4 = CSV.read("data_matlab/RndEx4.csv", DataFrame, header = 1)

M_NonStiff = CSV.read("data_matlab/BenchNonStiff.csv", DataFrame, header = 0)
M_Stiff = CSV.read("data_matlab/BenchStiff.csv", DataFrame, header = 0)
M_Harmonic = CSV.read("data_matlab/BenchHarmonic.csv", DataFrame, header = 0)
M_SIR = CSV.read("data_matlab/BenchSIR.csv", DataFrame, header = 0)
M_LV = CSV.read("data_matlab/BenchLV.csv", DataFrame, header = 0)

# Data from Julia
J_Nonstiff_E1 = CSV.read("data_Julia/NonStiff_E1.csv", DataFrame, header = 1)
J_Nonstiff_E2 = CSV.read("data_Julia/NonStiff_E2.csv", DataFrame, header = 1)
J_Nonstiff_T1 = CSV.read("data_Julia/NonStiff_T1.csv", DataFrame, header = 1)
J_Nonstiff_T2 = CSV.read("data_Julia/NonStiff_T2.csv", DataFrame, header = 1)

J_stiff_E1 = CSV.read("data_Julia/Stiff_E1.csv", DataFrame, header = 1)
J_stiff_E2 = CSV.read("data_Julia/Stiff_E2.csv", DataFrame, header = 1)
J_stiff_T1 = CSV.read("data_Julia/Stiff_T1.csv", DataFrame, header = 1)
J_stiff_T2 = CSV.read("data_Julia/Stiff_T2.csv", DataFrame, header = 1)

J_Harmonic_E1 = CSV.read("data_Julia/Harmonic_E1.csv", DataFrame, header = 1)
J_Harmonic_E2 = CSV.read("data_Julia/Harmonic_E2.csv", DataFrame, header = 1)
J_Harmonic_T1 = CSV.read("data_Julia/Harmonic_T1.csv", DataFrame, header = 1)
J_Harmonic_T2 = CSV.read("data_Julia/Harmonic_T2.csv", DataFrame, header = 1)

J_LV_E1 = CSV.read("data_Julia/LV_E1.csv", DataFrame, header = 1)
J_LV_E2 = CSV.read("data_Julia/LV_E2.csv", DataFrame, header = 1)
J_LV_T1 = CSV.read("data_Julia/LV_T1.csv", DataFrame, header = 1)
J_LV_T2 = CSV.read("data_Julia/LV_T2.csv", DataFrame, header = 1)
DynLV = CSV.read("data_Julia/DynLV.csv", DataFrame, header = 1)

J_SIR_E1 = CSV.read("data_Julia/SIR_E1.csv", DataFrame, header = 1)
J_SIR_E2 = CSV.read("data_Julia/SIR_E2.csv", DataFrame, header = 1)
J_SIR_T1 = CSV.read("data_Julia/SIR_T1.csv", DataFrame, header = 1)
J_SIR_T2 = CSV.read("data_Julia/SIR_T2.csv", DataFrame, header = 1)
# J_SIR_E3 = CSV.read("SIR_E3.csv", DataFrame, header = 1)
J_SIR_E4 = CSV.read("data_Julia/SIR_E4.csv", DataFrame, header = 1)
# J_SIR_T3 = CSV.read("SIR_T3.csv", DataFrame, header = 1)
J_SIR_T4 = CSV.read("data_Julia/SIR_T4.csv", DataFrame, header = 1)
J_SIR_E5 = CSV.read("data_Julia/SIR_E5.csv", DataFrame, header = 1)
J_SIR_E6 = CSV.read("data_Julia/SIR_E6.csv", DataFrame, header = 1)
J_SIR_T5 = CSV.read("data_Julia/SIR_T5.csv", DataFrame, header = 1)
J_SIR_T6 = CSV.read("data_Julia/SIR_T6.csv", DataFrame, header = 1)
J_SIR_E7 = CSV.read("data_Julia/SIR_E7.csv", DataFrame, header = 1)
J_SIR_E8 = CSV.read("data_Julia/SIR_E8.csv", DataFrame, header = 1)
J_SIR_E9 = CSV.read("data_Julia/SIR_E9.csv", DataFrame, header = 1)
J_SIR_T7 = CSV.read("data_Julia/SIR_T7.csv", DataFrame, header = 1)
J_SIR_T8 = CSV.read("data_Julia/SIR_T8.csv", DataFrame, header = 1)
J_SIR_T9 = CSV.read("data_Julia/SIR_T9.csv", DataFrame, header = 1)

J_tRnd1=CSV.read("data_Julia/tRnd1.csv", DataFrame, header = 1)
J_tRnd2=CSV.read("data_Julia/tRnd2.csv", DataFrame, header = 1)
J_ErrRnd1=CSV.read("data_Julia/ErrRnd1.csv", DataFrame, header = 1)
J_ErrRnd2=CSV.read("data_Julia/ErrRnd2.csv", DataFrame, header = 1)


## plot Examples

plot(J_Nonstiff_T1[:,1], J_Nonstiff_E1[:,1], xscale = :log, yscale = :log,
        legend_position=:bottomleft,
     label = "J1",c="firebrick3", shape = :circle, thickness_scaling = 1)
 plot!(J_Nonstiff_T2[:,1], J_Nonstiff_E2[:,1],label = "J2", c="hotpink",shape = :rect)
 plot!(M_NonStiff[:, 2], M_NonStiff[:, 6],label = "M1", shape = :circle, c="royalblue3")
 plot!(M_NonStiff[:, 4], M_NonStiff[:, 8],label = "M2", shape = :rect, c="skyblue2",
    title = "(a)", titleloc = :left, titlefont = font(10),legend=:false)
    plot!(M_NonStiff[:, 3], M_NonStiff[:, 7], label = "M3", shape = :diamond, c="cyan3")
     p1=plot!(M_NonStiff[:, 1], M_NonStiff[:, 5],label = "M4",shape = :rtriangle, c="mediumblue",
     ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",xlabel="Execution time (Sec)")

plot(J_stiff_T1[2:end,1], J_stiff_E1[2:end,1], xscale = :log, yscale = :log,
         legend_position=:bottomleft,
     label = "J1",c="firebrick3", shape = :circle, thickness_scaling = 1)
 plot!(J_stiff_T2[:,1], J_stiff_E2[:,1],label = "J2",c="hotpink", shape = :rect)
 plot!(M_Stiff[2:end, 2], M_Stiff[2:end, 6], label = "M1",c="royalblue3", shape = :circle)
 plot!(M_Stiff[:, 4], M_Stiff[:, 8], label = "M2",c="skyblue2", shape = :rect,legendfontsize=6,
 ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",xlabel="Execution time (Sec)",
 title = "(b)", titleloc = :left, legendtitle="Method",
 titlefont = font(10),legendposition=:outerbottom)
 plot!(M_Stiff[:, 3], M_Stiff[:, 7], label = "M3",c="cyan3", shape = :diamond)
 p2=plot!(M_Stiff[2:end, 1], M_Stiff[2:end, 5], label = "M4",c="mediumblue",shape = :rtriangle)

plot(J_Harmonic_T1[:,1], J_Harmonic_E1[:,1], xscale = :log, yscale = :log,
         legend_position=:bottomleft,
     label = "J1",c="firebrick3", shape = :circle, thickness_scaling = 1)
 plot!(J_Harmonic_T2[:,1], J_Harmonic_E2[:,1],label = "J2",c="hotpink", shape = :rect)
 plot!(M_Harmonic[:, 2], M_Harmonic[:, 6], label = "M1",c="royalblue3", shape = :circle)
 plot!(M_Harmonic[:, 4], M_Harmonic[:, 8], label = "M2",c="skyblue2", shape = :rect,
    ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",xlabel="Execution time (Sec)",
    title = "(c)", titleloc = :left, titlefont = font(10),legend=:false)
    plot!(M_Harmonic[:, 3], M_Harmonic[:, 7], label = "M3",c="cyan3", shape = :diamond)
    p3=plot!(M_Harmonic[:, 1], M_Harmonic[:, 5], label = "M4",c="mediumblue",shape = :rtriangle)

plot(J_SIR_T1[:,1], J_SIR_E1[:,1], xscale = :log, yscale = :log,
         legend_position=:bottomleft,
     label = "J1",c="firebrick3", shape = :circle, thickness_scaling = 1)
     plot!(J_SIR_T2[:,1], J_SIR_E2[:,1],label = "J2",c="hotpink", shape = :rect)
     plot!(J_SIR_T9[:,1], J_SIR_E9[:,1],   label = "J3", c="darkorange", shape = :circle)
     plot!(J_SIR_T4[:,1], J_SIR_E4[:,1], label = "J4", c="darkorange", shape = :hexagon)
     plot!(J_SIR_T5[:,1], J_SIR_E5[:,1], xscale = :log, yscale = :log,label = "J5", shape = :star5,c="darkorange")
     plot!(J_SIR_T6[:,1], J_SIR_E6[:,1],  label = "J6", c="darkorange", shape = :utriangle)
     plot!(J_SIR_T7[:,1], J_SIR_E7[:,1],  label = "J7", c="darkorange", shape = :dtriangle)
     plot!(J_SIR_T8[:,1], J_SIR_E8[:,1], label = "J8", c="darkorange", shape = :pentagon,
     ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",xlabel="Execution time (Sec)", legendfontsize=6,
     legendposition=:outerbottomright,title = "(a)", titleloc = :left, titlefont = font(10))
    plot!(M_SIR[:, 2], M_SIR[:, 6], label = "M1",c="royalblue3", shape = :circle)
    plot!(M_SIR[:, 4], M_SIR[:, 8], label = "M2",c="skyblue2", shape = :rect,legendtitle="Method")
    plot!(M_SIR[:, 3], M_SIR[:, 7], label = "M3",c="cyan3", shape = :diamond)
    p4=plot!(M_SIR[:, 1], M_SIR[:, 5], label = "M4",c="mediumblue",shape = :rtriangle)


plot(J_LV_T1[:,1], J_LV_E1[:,1], xscale = :log, yscale = :log,
         legend_position=:bottomleft,
     label = "J1",c="firebrick3", shape = :circle, thickness_scaling = 1)
     plot!(J_LV_T2[:,1], J_LV_E2[:,1],label = "J2",c="hotpink", shape = :rect)
     plot!(M_LV[:, 2], M_LV[:, 6], label = "M1",c="royalblue3", shape = :circle)
     plot!(M_LV[:, 4], M_LV[:, 8], label = "M2",c="skyblue2", shape = :rect,
     ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",xlabel="Execution time (Sec)",
     title = "(b)", titleloc = :left, titlefont = font(10)   ,legend=:false)
     plot!(M_LV[:, 3], M_LV[:, 7], label = "M3",c="cyan3", shape = :diamond)
     p5=plot!(M_LV[:, 1], M_LV[:, 5], label = "M4",c="mediumblue",shape = :rtriangle)

p6=plot(DynLV[:,1],Matrix(DynLV[:,2:4]),xlabel="Time", ylabel="Abundance of species" ,
    thickness_scaling = 1 , labels=["X1" "X2" "X3"],
    title = "(c)", titleloc = :left, titlefont = font(10))

# P=plot(p1, p2, layout = (2, 1),legend_position= (.6,.7) , size = (500, 500))
# P=plot(p1, p2, p3, layout = (2,2) , size = (1000, 1000))

# l = @layout [a b; c{.8w} _]
l = @layout [[grid(2,1)] b{.5w}]
# plot(p1, p2, p3,p4, layout = grid(3, 2, widths=[.14 ,0.4, 4,.4]) )
Plt1D=plot(p1, p3, p2, layout = l, size = (800, 600))


l = @layout [b{.6h}; grid(1,2)]
# plot(p1, p2, p3,p4, layout = grid(3, 2, widths=[.14 ,0.4, 4,.4]) )
PltMD=plot(p4, p5, p6, layout = l,size = (800, 600))


#### plot randoms

# scatter(Vector(J_tRnd1[:,1]), Vector(J_ErrRnd1[:,1]), xscale = :log, yscale = :log, linewidth = 3, markersize = 5,
#      label = "J1",c="firebrick3", shape = :circle, xlabel="Execution time (Sec)", ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",
#      thickness_scaling = 1,legend_position= :bottomleft, fc=:transparent,framestyle=:box)
#      scatter!(Vector(J_tRnd2[:,1]), Vector(J_ErrRnd2[:,1]),linewidth = 3, markersize = 5,label = "J2",c="firebrick3", shape = :rect)
#
#     Mt1=vcat(M_Ex1.Bench1_1,M_Ex2.Bench2_1,M_Ex3.Bench3_1,M_Ex4.Bench4_1)
#     Merr1=vcat(M_Ex1.Bench1_5,M_Ex2.Bench2_5,M_Ex3.Bench3_5,M_Ex4.Bench4_5)
#
#     Mt2=vcat(M_Ex1.Bench1_2,M_Ex2.Bench2_2,M_Ex3.Bench3_2,M_Ex4.Bench4_2)
#     Merr2=vcat(M_Ex1.Bench1_6,M_Ex2.Bench2_6,M_Ex3.Bench3_6,M_Ex4.Bench4_6)
#
#     Mt3=vcat(M_Ex1.Bench1_3,M_Ex2.Bench2_3,M_Ex3.Bench3_3,M_Ex4.Bench4_3)
#     Merr3=vcat(M_Ex1.Bench1_7,M_Ex2.Bench2_7,M_Ex3.Bench3_7,M_Ex4.Bench4_7)
#
#     Mt4=vcat(M_Ex1.Bench1_4,M_Ex2.Bench2_4,M_Ex3.Bench3_4,M_Ex4.Bench4_4)
#     Merr4=vcat(M_Ex1.Bench1_8,M_Ex2.Bench2_8,M_Ex3.Bench3_8,M_Ex4.Bench4_8)
#
#     scatter!(Mt1, Merr1, linewidth = 3, markersize = 5,label = "M4",c="royalblue2",shape = :rtriangle)
#     scatter!(Mt3, Merr3, linewidth = 3, markersize = 5,label = "M3",c="royalblue2", shape = :diamond)
#     scatter!(Mt2, Merr2, linewidth = 3, markersize = 5,label = "M1",c="royalblue2", shape = :circle)
#     p7=scatter!(Mt4, Merr4, linewidth = 3, markersize = 5,label = "M2",c="royalblue2", shape = :rect,
#     title = "(a)", titleloc = :left, titlefont = font(10))
#
# X = ["J1", "J2","M4","M3","M1","M2"]
#     Y=[Vector(J_ErrRnd1[:,1]),Vector(J_ErrRnd2[:,1]),Merr1,Merr2,Merr3,Merr4]
#     X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
#     df = DataFrame(X = X2, Y = Y)
#     p8 = @df df boxplot(:X, :Y, c=:black, fillcolor=:white, legend=false,
#     title = "(b)", titleloc = :left, titlefont = font(10),
#     yscale=:log, ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2")
#
# X = ["J1", "J2","M4","M3","M1","M2"]
#         Y=[Vector(J_tRnd1[:,1]),Vector(J_tRnd2[:,1]),Mt1,Mt2,Mt3,Mt4]
#         X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
#         df = DataFrame(X = X2, Y = Y)
#         p9 = @df df boxplot(:X, :Y, c=:black, fillcolor=:white, legend=false,
#         title = "(c)", titleloc = :left, titlefont = font(10),
#         yscale=:log, ylabel="Execution time (Sec)")
#
# l = @layout [b{.5w} grid(2,1)]
# PltRnd=plot(p7, p8, p9, layout = l,size = (800, 500))



##

Mt1=vcat(M_Ex1.Bench1_1,M_Ex2.Bench2_1,M_Ex3.Bench3_1,M_Ex4.Bench4_1)
    Merr1=vcat(M_Ex1.Bench1_5,M_Ex2.Bench2_5,M_Ex3.Bench3_5,M_Ex4.Bench4_5)

    Mt2=vcat(M_Ex1.Bench1_2,M_Ex2.Bench2_2,M_Ex3.Bench3_2,M_Ex4.Bench4_2)
    Merr2=vcat(M_Ex1.Bench1_6,M_Ex2.Bench2_6,M_Ex3.Bench3_6,M_Ex4.Bench4_6)

    Mt3=vcat(M_Ex1.Bench1_3,M_Ex2.Bench2_3,M_Ex3.Bench3_3,M_Ex4.Bench4_3)
    Merr3=vcat(M_Ex1.Bench1_7,M_Ex2.Bench2_7,M_Ex3.Bench3_7,M_Ex4.Bench4_7)

    Mt4=vcat(M_Ex1.Bench1_4,M_Ex2.Bench2_4,M_Ex3.Bench3_4,M_Ex4.Bench4_4)
    Merr4=vcat(M_Ex1.Bench1_8,M_Ex2.Bench2_8,M_Ex3.Bench3_8,M_Ex4.Bench4_8)

    #plot linear lines
# td=log.(J_tRnd1[:,1])
#     ed=log.(J_ErrRnd1[:,1])
#     data = DataFrame(x=td, y=ed);
#
#     model = lm(@formula(y ~ x), data)
#     pred = DataFrame(x = minimum(td):maximum(td))
#     pr = predict(model, pred, interval = :confidence, level = 0.0)
#
#     plot(;xscale = :log, yscale = :log)
#     plot!(exp.(pred.x), exp.(pr.prediction),  linewidth=3, c="firebrick3",
#             ribbon = (exp.(pr.prediction) .- exp.(pr.lower), exp.(pr.upper) .- exp.(pr.prediction)))
#             td=log.(J_tRnd2[:,1])
#     ed=log.(J_ErrRnd2[:,1])
#     data = DataFrame(x=td, y=ed);
#
#     model = lm(@formula(y ~ x), data)
#     pred = DataFrame(x = minimum(td):maximum(td))
#     pr = predict(model, pred, interval = :confidence, level = 0.0)
#     plot!(exp.(pred.x), exp.(pr.prediction),  linewidth=3, c="hotpink",
#             ribbon = (exp.(pr.prediction) .- exp.(pr.lower), exp.(pr.upper) .- exp.(pr.prediction)))
#             td=log.(Mt2)
#     ed=log.(Merr2)
#     data = DataFrame(x=td, y=ed);
#
#     model = lm(@formula(y ~ x), data)
#     pred = DataFrame(x = minimum(td):maximum(td))
#     pr = predict(model, pred, interval = :confidence, level = 0.0)
#     plot!(exp.(pred.x), exp.(pr.prediction),  linewidth=3, c="royalblue3",
#             ribbon = (exp.(pr.prediction) .- exp.(pr.lower), exp.(pr.upper) .- exp.(pr.prediction)))
#             td=log.(Mt4)
#     ed=log.(Merr4)
#     data = DataFrame(x=td, y=ed);
#
#     model = lm(@formula(y ~ x), data)
#     pred = DataFrame(x = minimum(td):maximum(td))
#     pr = predict(model, pred, interval = :confidence, level = 0.0)
#     plot!(exp.(pred.x), exp.(pr.prediction),  linewidth=3, c="skyblue2",
#             ribbon = (exp.(pr.prediction) .- exp.(pr.lower), exp.(pr.upper) .- exp.(pr.prediction)))
#             td=log.(Mt3)
#     ed=log.(Merr3)
#     data = DataFrame(x=td, y=ed);
#
#     model = lm(@formula(y ~ x), data)
#     pred = DataFrame(x = minimum(td):maximum(td))
#     pr = predict(model, pred, interval = :confidence, level = 0.0)
#     plot!(exp.(pred.x), exp.(pr.prediction),  linewidth=3, c="cyan3",
#             ribbon = (exp.(pr.prediction) .- exp.(pr.lower), exp.(pr.upper) .- exp.(pr.prediction)))
#             td=log.(Mt1)
#     ed=log.(Merr1)
#     data = DataFrame(x=td, y=ed);
#
#     model = lm(@formula(y ~ x), data)
#     pred = DataFrame(x = minimum(td):maximum(td))
#     pr = predict(model, pred, interval = :confidence, level = 0.0)
#     plot!(exp.(pred.x), exp.(pr.prediction),  linewidth=3, c="mediumblue",
#             ribbon = (exp.(pr.prediction) .- exp.(pr.lower), exp.(pr.upper) .- exp.(pr.prediction)))
Y=fitlinear(log.(J_tRnd1[:,1]),log.(J_ErrRnd1[:,1]))
    plot(exp.(Y.x),exp.(Y.y),c="firebrick3",linewidth = 2,label=:false)

    Y=fitlinear(log.(J_tRnd2[:,1]),log.(J_ErrRnd2[:,1]))
    plot!(exp.(Y.x),exp.(Y.y),color = :hotpink,linewidth = 2,label=:false)

    Y=fitlinear(log.(Mt2),log.(Merr2))
    plot!(exp.(Y.x),exp.(Y.y),color = "royalblue3",linewidth = 2,label=:false)

    Y=fitlinear(log.(Mt4),log.(Merr4))
    plot!(exp.(Y.x),exp.(Y.y), color = :skyblue2,linewidth = 2,label=:false)

    Y=fitlinear(log.(Mt3),log.(Merr3))
    plot!(exp.(Y.x),exp.(Y.y), color = :cyan3,linewidth = 2,label=:false)

    Y=fitlinear(log.(Mt1),log.(Merr1))
    plot!(exp.(Y.x),exp.(Y.y), color = :mediumblue,linewidth = 2,label=:false)

scatter!(Vector(J_tRnd1[:,1]), Vector(J_ErrRnd1[:,1]), markerstrokewidth=0,xscale = :log, yscale = :log,
         label = "J1",c="firebrick3", shape = :circle, xlabel="Execution time (Sec)", ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2",
         thickness_scaling = 1,legend_position= :bottomleft, fc=:transparent)
    scatter!(Vector(J_tRnd2[:,1]), Vector(J_ErrRnd2[:,1]),markerstrokewidth=0, label = "J2",c="hotpink", shape = :rect)

    scatter!(Mt2, Merr2, label = "M1",c="royalblue3", shape = :circle,markerstrokewidth=0)
    scatter!(Mt4, Merr4,  label = "M2",c="skyblue2", shape = :rect,markerstrokewidth=0,
    titleloc = :left, titlefont = font(10), legendtitle="Method", legendtitlefont = font(9))
    scatter!(Mt3, Merr3, label = "M3",c="cyan3", shape = :diamond,markerstrokewidth=0)
    p10=scatter!(Mt1, Merr1,  label = "M4",c="mediumblue",shape = :rtriangle,markerstrokewidth=0, ms=6,
    title = "(a)", titleloc = :left, titlefont = font(10))

ColorSet=[:firebrick3 :hotpink :royalblue3 :skyblue2 :cyan3 :mediumblue]
    X = ["J1", "J2","M1","M2","M3","M4"]
    Y=[Vector(J_ErrRnd1[:,1]),Vector(J_ErrRnd2[:,1]),Merr2,Merr4,Merr3,Merr1]
    X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
    df = DataFrame(X = X2, Y = Y)
    p8 = @df df boxplot(:X, :Y , c=ColorSet,fillcolor=ColorSet, legend=false,markerstrokewidth=0,
    title = "(b)", titleloc = :left, titlefont = font(10),
    yscale=:log, ylabel=L"\textrm{Error}= \Vert x - \overline{x} \;\Vert _2")

X = ["J1", "J2","M1","M2","M3","M4"]
        Y=[Vector(J_tRnd1[:,1]),Vector(J_tRnd2[:,1]),Mt2,Mt4,Mt3,Mt1]
        X2 = [fill(x,length(y)) for (x,y) in zip(X,Y)]
        df = DataFrame(X = X2, Y = Y)
        p9 = @df df boxplot(:X, :Y, c=ColorSet,fillcolor=ColorSet, legend=false,markerstrokewidth=0,
        title = "(c)", titleloc = :left, titlefont = font(10),
        yscale=:log, ylabel="Execution time (Sec)")


l = @layout [b{.5w} grid(2,1)]
PltRnd=plot(p10, p8, p9, layout = l,size = (800, 500))


savefig(Plt1D,"Plt1D.svg")
savefig(PltMD,"PltMD.svg")
savefig(PltRnd,"PltRnd.svg")
savefig(Plt1D,"Plt1D.png")
savefig(PltMD,"PltMD.png")
savefig(PltRnd,"PltRnd.png")
