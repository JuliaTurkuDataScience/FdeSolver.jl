push!(LOAD_PATH, "../src/")
ENV["GKS_WSTYPE"]=100
using FdeSolver, Documenter

makedocs(
         format=Documenter.HTML(;
         canonical="https://github.com/JuliaTurkuDataScience/FdeSolver.jl"
         ),
         sitename = "FdeSolver.jl",
         modules  = [FdeSolver],
         pages=[
                "Home" => "intro.md"
                "Manual" => "index.md"
                "Examples" => "examples.md"
               ])

deploydocs(;
        repo="github.com/JuliaTurkuDataScience/FdeSolver.jl",
)
