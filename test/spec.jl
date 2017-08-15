@testset "Spec" begin

  uc = newunitcell([1.0 0.0; 0.0 1.0])
  addorbital!(uc, "A", FractCoord([0,0], [0.5, 0.0]))
  addorbital!(uc, "B", FractCoord([0,0], [0.0, 0.5]))

  let
    hopping = hoppingbycarte(uc, 1.0, "A", "B", [0.5, 0.0], [0.0, 0.5])
    @test islocal(hopping)
    @test islocal(localize(hopping))
    @test hopping ≂ localize(hopping)
  end
  let
    hopping = hoppingbycarte(uc, 1.0, "A", "B", [1.5, 0.0], [1.0, 0.5])
    @test !islocal(hopping)
    @test islocal(localize(hopping))
    @test !(hopping ≂ localize(hopping))
  end
end
