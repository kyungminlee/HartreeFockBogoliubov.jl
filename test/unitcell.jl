@testset "UnitCell" begin
  latticevectors = [0.5 0.0; 0.0 1.0]

  @testset "Constructor" begin
    uc = newunitcell(latticevectors)
    @test isapprox(uc.latticevectors, latticevectors)
    @test isapprox(uc.reducedreciprocallatticevectors, [2.0 0.0; 0.0 1.0])
    @test isapprox(uc.reciprocallatticevectors, [4*pi 0.0; 0.0 2*pi])
  end

  @testset "Methods" begin
    uc = newunitcell(latticevectors)
    fc1 = FractCoord([0, 0], [0.5, 0.0])
    fc2 = FractCoord([0, 0], [0.0, 0.5])
    addorbital!(uc, "Ox", fc1)
    addorbital!(uc, "Oy", fc2)
    @test numorbital(uc) == 2
    @test hasorbital(uc, "Ox")
    @test hasorbital(uc, "Oy")
    @test !hasorbital(uc, "Oz")
    @test getorbitalindex(uc, "Ox") == 1
    @test getorbitalindex(uc, "Oy") == 2
    @test getorbital(uc, "Ox") == ("Ox", fc1)
    @test getorbital(uc, "Oy") == ("Oy", fc2)
    @test getorbitalcoord(uc, "Ox") == fc1
    @test getorbitalcoord(uc, "Oy") == fc2
    @test getorbitalindexcoord(uc, "Ox") == (1, fc1)
    @test getorbitalindexcoord(uc, "Oy") == (2, fc2)

    @test getorbital(uc, 1) == ("Ox", fc1)
    @test getorbital(uc, 2) == ("Oy", fc2)
    @test getorbitalname(uc, 1) == "Ox"
    @test getorbitalname(uc, 2) == "Oy"
    @test getorbitalcoord(uc, 1) == fc1
    @test getorbitalcoord(uc, 2) == fc2
  end

  @testset "Type" begin
    uc = newunitcell(latticevectors; OrbitalType=typeof((:up, "A")))
    fc = Dict("Ox" => FractCoord([0, 0], [0.5, 0.0]),
              "Oy" => FractCoord([0, 0], [0.0, 0.5]))
    for orbital in ["Ox", "Oy"]
      for spin in [:up, :dn]
        addorbital!(uc, (spin, orbital), fc[orbital])
      end
    end
    @test getorbitalname(uc, 1) == (:up, "Ox")
    @test getorbitalname(uc, 2) == (:dn, "Ox")

    @show uc

  end
end
