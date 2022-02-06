using ReactionEngine
using Documenter

DocMeta.setdocmeta!(ReactionEngine, :DocTestSetup, :(using ReactionEngine); recursive=true)

makedocs(;
    modules=[ReactionEngine],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/ReactionEngine.jl/blob/{commit}{path}#{line}",
    sitename="ReactionEngine.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/ReactionEngine.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/ReactionEngine.jl",
    devbranch="main",
)
