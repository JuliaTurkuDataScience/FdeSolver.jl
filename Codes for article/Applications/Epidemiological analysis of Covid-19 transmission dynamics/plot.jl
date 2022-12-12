#plotting
using Plots,StatsPlots, StatsPlots.PlotMeasures

P=plot(plSpain,plPortugal, layout = grid(2,1), size=(700,650))

savefig(P,"CovidSpnPrt.svg")
