using Gm_ID_kit
using Documenter

DocMeta.setdocmeta!(Gm_ID_kit, :DocTestSetup, :(using Gm_ID_kit); recursive=true)

makedocs(;
    modules=[Gm_ID_kit],
    authors="islam omar ahmed <io1131@fayoum.edu.eg> and contributors",
    repo="https://github.com/Islam0mar/Gm_ID_kit.jl/blob/{commit}{path}#{line}",
    sitename="Gm_ID_kit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Islam0mar.github.io/Gm_ID_kit.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Islam0mar/Gm_ID_kit.jl",
    devbranch="main",
)
