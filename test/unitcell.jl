@testset "UnitCell" begin
  latticevectors = [0.5 0.0; 0.0 1.0]

  @testset "Constructor" begin
    uc = make_unitcell(latticevectors)
    make_unitcell([[1.0,0.0], [0.0,1.0]])
    @test isapprox(uc.latticevectors, latticevectors)
    @test isapprox(uc.reducedreciprocallatticevectors, [2.0 0.0; 0.0 1.0])
    @test isapprox(uc.reciprocallatticevectors, [4*pi 0.0; 0.0 2*pi])
  end

  @testset "Constructor Exceptions" begin
    @test_throws ArgumentError make_unitcell(latticevectors; OrbitalType=Int)
    @test_throws ArgumentError make_unitcell([1.0 0.0;]; OrbitalType=String)
    @test_throws ArgumentError make_unitcell([1.0 0.0; 1.0 0.0]; OrbitalType=String)
    @test_throws ArgumentError make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String, tol=-1E-8)
  end

  @testset "Methods" begin
    uc = make_unitcell(latticevectors)
    fc1 = FractCoord([0, 0], [0.5, 0.0])
    fc2 = FractCoord([0, 0], [0.0, 0.5])
    index1 = addorbital!(uc, "Ox", fc1)
    index2 = addorbital!(uc, "Oy", fc2)
    @test index1 == 1
    @test index2 == 2

    @test dimension(uc) == 2
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

  @testset "Methods Exceptions" begin
    uc = make_unitcell(latticevectors; OrbitalType=String)
    addorbital!(uc, "Ox", FractCoord([0, 0], [0.0, 0.0]))
    @test_throws ArgumentError addorbital!(uc, "Oy", FractCoord([0], [0.5]))
    @test_throws ArgumentError addorbital!(uc, "Ox", FractCoord([0, 0], [0.0, 0.0]))
  end

  @testset "Type" begin
    uc = make_unitcell(latticevectors; OrbitalType=typeof((:up, "A")))
    fc = Dict("Ox" => FractCoord([0, 0], [0.5, 0.0]),
              "Oy" => FractCoord([0, 0], [0.0, 0.5]))
    for orbital in ["Ox", "Oy"]
      for spin in [:up, :dn]
        addorbital!(uc, (spin, orbital), fc[orbital])
      end
    end
    @test getorbitalname(uc, 1) == (:up, "Ox")
    @test getorbitalname(uc, 2) == (:dn, "Ox")
  end

  @testset "Conversion fract2carte/carte2fract" begin
    latticevectors = [0.5 0.0; 0.0 1.0]
    uc = make_unitcell(latticevectors)

    rawfractcoord = [-1.2, 1.5]

    correctfractcoord = FractCoord([-2, 1], [0.8, 0.5])
    correctcartecoord = [-0.6, 1.5]

    fractcoord = FractCoord(rawfractcoord)
    cartecoord = fract2carte(uc, fractcoord)
    newfractcoord = carte2fract(uc, cartecoord)
    @test isapprox(fractcoord, correctfractcoord)
    @test isapprox(cartecoord, correctcartecoord)
    @test isapprox(newfractcoord, correctfractcoord)
  end




  @testset "momentumgrid" begin
    uc = make_unitcell([1.0 0.0; 0.0 1.0]; OrbitalType=String)
    kg = momentumgrid(uc, [2,3])
    for i in 1:2
        for j in 1:3
            @test isapprox(kg[i,j], [2 * pi * (i-1) / 2, 2 * pi * (j-1) / 3])
        end
    end
  end
end
