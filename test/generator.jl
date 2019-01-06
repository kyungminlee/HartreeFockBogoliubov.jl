using LinearAlgebra

#=
@testset "Generator" begin
  uc = UnitCell([1.0 0.0; 0.0 1.0])
  addorbital!(uc, "A", FractCoord([0,0], [0.5, 0.0]))
  addorbital!(uc, "B", FractCoord([0,0], [0.0, 0.5]))

  hopspec_dia = Spec.HoppingDiagonal(1.0, "A", [0, 0])
  hopspec_off = Spec.HoppingOffdiagonal(1.0 + 0.3im, "A", "B", [0, 0], [0, 0])

  hopembed_dia = Embed.EmbedHoppingDiagonal(uc, hopspec_dia)
  hopembed_off = Embed.EmbedHoppingOffdiagonal(uc, hopspec_off)

  func = Generator.hopping_inplace(uc, hopembed_dia)
  out = zeros(ComplexF64, (2,2))

  @show out
  func([0.0, 0.0], out)
  @show out

  func = Generator.hopping_inplace(uc, hopembed_off)
  out = zeros(ComplexF64, (2,2))
  @show out
  func([0.0, 0.0], out)
  @show out

  out = zeros(ComplexF64, (2,2))
  @show out
  func([0.1, 0.0], out)
  @show out
end
=#

#=
 Dispersion(1 Dimension, 1 Orbital)
=#


@testset "disp_1D1O" begin
  uc = make_unitcell(1.0)
  addorbital!(uc, "A", FractCoord([0], [0.0]))
  t0 = 100.0
  t1 = 10.0 + 0.1im
  t2 = 1.0 + 0.01im

  hamspec = FullHamiltonian(uc)
  addhopping!(hamspec, hoppingbycarte(uc, t0, "A", [0.0]))
  addhopping!(hamspec, hoppingbycarte(uc, t1, "A", "A", [0.0], [1.0]))
  addhopping!(hamspec, hoppingbycarte(uc, t2, "A", "A", [0.0], [2.0]))

  func1 = Generator.hopping_inplace(uc, hamspec.hoppings_diagonal, hamspec.hoppings_offdiagonal)

  for kx in range(0.0, stop=2.0*pi, length=16+1)
    out1 = zeros(ComplexF64, (1,1))
    out2 = zeros(ComplexF64, (1,1))
    func1([kx], out1)
    ek = (t0
          + t1 * cis(kx) + conj(t1) * cis(-kx)
          + t2 * cis(2*kx) + conj(t2) * cis(-2*kx) )
    #@show out[1,1]
    #@show ek
    @test isapprox(out1[1,1], ek)
  end
end

#=
    Dispersion (2D, Graphene)
=#
@testset "graphene" begin
  a1 = [ sqrt(3.0) * 0.5, 0.5]
  a2 = [-sqrt(3.0) * 0.5, 0.5]
  a3 = [ 0.0, -1.0]
  b1 = a1 - a2
  b2 = a2 - a3
  b3 = a3 - a1
  uc = make_unitcell([b1 b3])
  addorbital!(uc, "A", carte2fract(uc, [0.0, 0.0]))
  addorbital!(uc, "B", carte2fract(uc, [0.0,-1.0]))
  #@show uc

  hopspecs = Spec.HoppingOffdiagonal[]
  t1 = 1.0 + 0.1im
  for r in [a1, a2, a3]
    push!( hopspecs, Spec.hoppingbycarte(uc, t1, "A", "B", [0.0, 0.0], r) )
  end

  func1 = Generator.hopping_inplace(uc, Spec.HoppingDiagonal[], hopspecs)

  for kx in range(4.0, stop=4.0, length=8)
    for ky in range(-4.0, stop=4.0, length=8)
      k = [kx, ky]
      mat = begin
        phase = (cis(dot(k, a1)) + cis(dot(k, a2)) + cis(dot(k, a3)) )
        T1 = t1 * phase
        T1c = conj(T1)
        [  0  T1 ;  T1c  0 ]
      end

      out1 = zeros(ComplexF64, (2, 2))
      func1(k, out1)

      @test isapprox(out1, mat)
    end
  end
end


#=

@testset "Generator" begin

  uc = UnitCell([1.0 0.0; 0.0 1.0])
  addorbital!(uc, "A", FractCoord([0,0], [0.5, 0.0]))
  addorbital!(uc, "B", FractCoord([0,0], [0.0, 0.5]))

  hopspec = Spec.Hopping[
    Spec.HoppingDiagonal( 1.0, "A", [0, 0]),
    Spec.HoppingDiagonal(-1.0, "B", [0, 0]),
    Spec.HoppingOffdiagonal(1.0 + 0.1im, "A", "B", [0, 0], [0, 0]),
    #Spec.HoppingOffdiagonal(1.0 + 0.2im, "A", "B", [0, 0], [0, 1]),
    #Spec.HoppingOffdiagonal(1.0 + 0.3im, "A", "B", [0,-1], [0, 0]),
    Spec.HoppingOffdiagonal(1.0 + 0.4im, "A", "B", [0, 0], [-1,0]),
  ]
  hopembed = Embed.EmbedHopping[Embed.embed(uc, hop) for hop in hopspec]

  func = Generator.generate(uc, hopembed)
  out = zeros(ComplexF64, (2,2))
  @show out
  func([0.0, 0.0], out)
  @show out

  out = zeros(ComplexF64, (2,2))
  @show out
  func([0.0001, 0.0002], out)
  @show out

  out = zeros(ComplexF64, (2,2))
  @show out
  func([0.0001, 0.0002], out)
  @show out

end

=#
