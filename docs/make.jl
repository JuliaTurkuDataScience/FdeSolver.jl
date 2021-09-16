using Documenter, FdeSolver

DocMeta.setdocmeta!(FdeSolver, :DocTestSetup, :(using FdeSolver); recursive=true)

makedocs(;
    modules=[FdeSolver],
    authors="Moein Khalighi, Giulio Benedetti",
    repo="https://github.com/JuliaTurkuDataScience/FdeSolver.jl/blob/{commit}{path}#{line}",
    sitename="FdeSolver.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/JuliaTurkuDataScience/FdeSolver.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaTurkuDataScience/FdeSolver.jl",
)
