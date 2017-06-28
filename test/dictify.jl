using HartreeFockBogoliubov
using HartreeFockBogoliubov: Spec, Generator, HFB, Dictify
#using JSON
using DataStructures

function newhubbard(μ::Real, t::Real, U::Real, V::Real)
  unitcell = newunitcell(eye(2), OrbitalType=Symbol)
  addorbital!(unitcell, :UP, FractCoord([0, 0], [0.0, 0.0]))
  addorbital!(unitcell, :DN, FractCoord([0, 0], [0.0, 0.0]))

  hamspec = Spec.FullHamiltonian(unitcell)
  for sp in [:UP, :DN]
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -μ, sp, [0.0, 0.0]))
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, sp, sp, [0.0, 0.0], [ 1.0, 0.0]))
    Spec.addhopping!(hamspec, Spec.hoppingbycarte(unitcell, -t, sp, sp, [0.0, 0.0], [ 0.0, 1.0]))
  end

  Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, U, :UP, :DN, [0.0, 0.0], [0.0, 0.0]))

  for sp1 in [:UP, :DN], sp2 in [:UP, :DN]
    Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, sp1, sp2, [0.0, 0.0], [1.0, 0.0]))
    Spec.addinteraction!(hamspec, Spec.interactionbycarte(unitcell, V, sp1, sp2, [0.0, 0.0], [0.0, 1.0]))
  end

  return hamspec
end

using YAML
using JSON

if true
  o1 = newunitcell([1 2; 3 4.0])
  @show o1

  d1 = dictify(o1)
  @show d1

  j1 = JSON.json(d1)
  @show j1

  d2 = JSON.parse(j1; dicttype=OrderedDict)
  @show d2

  o2 = objectify(d2)
  @show o2

  exit()
end

begin
  μ = 0.1
  t = 0.2
  U = 0.3
  V = 0.4
  nx, ny = 4,4
  temperature = 0.0
  ham = newhubbard(μ,t,U,V)
  solver = HFBSolver(ham, [nx, ny], temperature)
  currentsolution = newhfbsolution(solver.hfbcomputer)

  #=
  dict = dictify(hamspec)
  obj = objectify(dict)
  @show dict
  @show obj
  =#

  #=
  foo = dictify(hamspec.unitcell.orbitals[1])
  @show foo
  bar = objectify(foo)
  @show bar
  =#

  println("Spec.Hamiltonian")
  println("================")
  dict = dictify(ham)
  @show dict

  obj = objectify(dict)
  @show obj

  println()
  println()

  println("HFBComputer")
  println("=================")
  dict = dictify(solver.hfbcomputer)
  @show dict

  obj = objectify(dict)
  @show obj

  println()
  println()

  println("HFBSolver")
  println("=================")
  dict = dictify(solver)
  @show JSON.json(dict)

  obj = objectify(dict)
  @show obj

  jj = JSON.json(dict)
  dict2 = JSON.parse(jj; dicttype=OrderedDict)

  obj2 = objectify(dict2)

  @show obj2

  objectify(false)

  open("solver.json", "w") do file
    dict = dictify(solver)
    write(file, JSON.json(dict))
  end
  open("computer.json", "w") do file
    dict = dictify(solver.hfbcomputer)
    write(file, JSON.json(dict))
  end


end
