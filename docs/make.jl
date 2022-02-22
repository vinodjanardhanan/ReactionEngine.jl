# using Pkg
# pkg"activate .."
using ReactionEngine
using Documenter
include("../src/Utils.jl")
include("../src/IdealGas.jl")
include("../src/Inspect.jl")
include("../src/Cstr.jl")
include("../src/Plug.jl")
include("../src/ThermoProbe.jl")
include("../src/Batch.jl")
include("../src/Equil.jl")
include("../src/Transport.jl")

DocMeta.setdocmeta!(ReactionEngine, :DocTestSetup, :(using ReactionEngine); recursive=true)

makedocs(;
    modules=[ReactionEngine, IdealGas, Inspect, Cstr, Plug, ThermoProbe, Batch, Equil, Transport],
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
        "Models"=>[
            "Inspect" => "Models/inspect.md"
            "Thermoprobe"=>"Models/tprobe.md"
            "Equilibrate"=>"Models/equil.md"
            "Batch"=>"Models/batch.md"
            "Cstr"=>"Models/cstr.md"
            "Plug"=>"Models/plug.md"
            "Sensitivity"=>"Models/sens.md"
        ],
        "Library"=>[
            "Thermodynamic Properties"=>"Library/thermo.md"
            "Surface Chemistry"=>"Library/schem.md"
            "Transport Properties"=>"Library/trans.md"
        ]

    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/ReactionEngine.jl",
    devbranch="main",
)
