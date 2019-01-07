using Documenter
using HartreeFockBogoliubov

makedocs(
	modules=[HartreeFockBogoliubov],
	doctest=false,
	clean=true,
	checkdocs=:all,
	format=Documenter.HTML(prettyurls=!("local" in ARGS)),
	authors="Kyungmin Lee",
	sitename="HartreeFockBogoliubov.jl",
	pages=[
			"Home" => "index.md",
			"Basic Theory" => [
				"Hartree-Fock-Bogoliubov Theory" => "basics/hartreefockbogoliubov.md",
				"Hartree-Fock-Bogoliubov Equations" => "basics/hfbequation.md",
				"Hartree-Fock-Bogoliubov Decomposition" => "basics/hfbdecomposition.md",
				"Momentum Space Formulation" => "basics/momentumspace.md",
				"Topological Invariants" => "basics/topology.md",
			],
			"User Guide" => [
				"Lattice" => "guide/lattice.md",
				"Hamiltonian" => "guide/hamiltonian.md",
				"HFB" => "guide/hfb.md",
				"Linearized Gap Equation" => "guide/linearizedgap.md",
				"Topology" => "guide/topology.md",
				"Example" => "guide/example.md",
			],
			hide("Internals" => "internals.md", 
				["internals/hartreefockbogoliubov.md",
			])
		]
	)

deploydocs(
  deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
  repo = "github.com/kyungminlee/HartreeFockBogoliubov.jl.git",
)
