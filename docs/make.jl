using Documenter
using HartreeFockBogoliubov

makedocs(sitename="HartreeFockBogoliubov.jl")

deploydocs(
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
  repo = "github.com/kyungminlee/HartreeFockBogoliubov.jl.git",
  julia = "1.0",
  osname = "linux",
)
