@testset "Spec" begin

  uc = newunitcell([1.0 0.0; 0.0 1.0])
  addorbital!(uc, "A", FractCoord([0,0], [0.5, 0.0]))
  addorbital!(uc, "B", FractCoord([0,0], [0.0, 0.5]))

  #hopspec_dia = SpecHoppingDiagonal(1.0, "A", [0, 0])
  #hopspec_off = SpecHoppingOffdiagonal(1.0, "A", "B", [0, 0], [0, 0])

  #@show hopspec_dia
  #@show hopspec_off

  #hopembed_dia = Embed.EmbedHoppingDiagonal(uc, hopspec_dia)
  #hopembed_off = Embed.EmbedHoppingOffdiagonal(uc, hopspec_off)

  #@show hopembed_dia
  #@show hopembed_off
end
