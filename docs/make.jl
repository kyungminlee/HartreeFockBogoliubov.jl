using Documenter
using HartreeFockBogoliubov

makedocs(
    authors="Kyungmin Lee",
    sitename="HartreeFockBogoliubov.jl",
    pages = Any[
        "HOME" => "index.md",
        "Basics" => Any[
            "basics/hartreefockbogoliubov.md",
            "basics/momentumspace.md"
        ]
    ]
)

deploydocs(
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
  repo = "github.com/kyungminlee/HartreeFockBogoliubov.jl.git",
  julia = "0.5",
)
