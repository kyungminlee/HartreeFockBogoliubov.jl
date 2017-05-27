push!(LOAD_PATH, "../src")

using Documenter
using HartreeFockBogoliubov


makedocs(
)


deploydocs(
  deps = Deps.pip("pygments", "mkdics", "python-markdown-math"),
  julia = "0.5"
)