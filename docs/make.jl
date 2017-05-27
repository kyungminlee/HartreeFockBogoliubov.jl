push!(LOAD_PATH, "../src")

using Documenter
using HartreeFockBogoliubov


makedocs(
)


deploydocs(
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
  repo = "https://github.com/kyungminlee/HartreeFockBogoliubov.jl.git",
  julia = "0.5",
)