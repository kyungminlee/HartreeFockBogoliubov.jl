using Documenter
using HartreeFockBogoliubov

makedocs()

deploydocs(
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
  repo = "github.com/kyungminlee/HartreeFockBogoliubov.jl.git",
  julia = "0.6",
  osname = "linux",
)
