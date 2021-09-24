push!(LOAD_PATH, "../src/")
using FdeSolver, Documenter

makedocs(
         sitename = "FdeSolver.jl",
         modules  = [FdeSolver],
         pages=[
                "Home" => "index.md"
               ])

deploydocs(;
        repo="github.com/JuliaTurkuDataScience/FdeSolver.jl",
)
