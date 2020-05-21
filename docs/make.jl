using Documenter
using Hodge

# Set the right metadata for doctests
DocMeta.setdocmeta!(Hodge, :DocTestSetup, :(using Hodge); recursive=true)

makedocs(
    sitename = "Hodge.jl",
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical  = "https://iagoleal.github.io/Hodge.jl/dev/",
        assets     = []
    ),
    modules = [Hodge],
    pages    = [
        "Introduction" => "index.md",
        "Reference" => [
            "Public" => "refs-api.md",
            "Internal" => "refs-simplextree.md"
        ]
    ]
)

# Deploy site to Github Pages
if !("local" in ARGS)
    deploydocs(
        repo = "github.com/iagoleal/Hodge.jl.git"
    )
end
