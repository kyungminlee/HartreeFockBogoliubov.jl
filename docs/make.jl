using Documenter
using HartreeFockBogoliubov

makedocs()

deploydocs(
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math", "mkdocs-cinder"),
  repo = "github.com/kyungminlee/HartreeFockBogoliubov.jl.git",
  julia = "0.6",
)
