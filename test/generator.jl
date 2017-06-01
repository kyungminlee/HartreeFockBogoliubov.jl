#=
@testset "Generator" begin
  uc = UnitCell([1.0 0.0; 0.0 1.0])
  addorbital!(uc, "A", FractCoord([0,0], [0.5, 0.0]))
  addorbital!(uc, "B", FractCoord([0,0], [0.0, 0.5]))

  hopspec_dia = Spec.HoppingDiagonal(1.0, "A", [0, 0])
  hopspec_off = Spec.HoppingOffdiagonal(1.0 + 0.3im, "A", "B", [0, 0], [0, 0])

  hopembed_dia = Embed.HoppingDiagonal(uc, hopspec_dia)
  hopembed_off = Embed.HoppingOffdiagonal(uc, hopspec_off)

  func = Generator.generatefast(uc, hopembed_dia)
  out = zeros(Complex128, (2,2))

  @show out
  func([0.0, 0.0], out)
  @show out

  func = Generator.generatefast(uc, hopembed_off)
  out = zeros(Complex128, (2,2))
  @show out
  func([0.0, 0.0], out)
  @show out

  out = zeros(Complex128, (2,2))
  @show out
  func([0.1, 0.0], out)
  @show out
end
=#

#=
 Dispersion(1 Dimension, 1 Orbital)
=#
@testset "disp_1D1O" begin
  uc = UnitCell([1.0])
  addorbital!(uc, "A", FractCoord([0], [0.0]))
  t0 = 100.0
  t1 = 10.0 + 0.1im
  t2 = 1.0 + 0.01im
  hopspec = Spec.Hopping[
    Spec.HoppingDiagonal( t0, "A", [0]),
    Spec.HoppingOffdiagonal( t1, "A", "A", [0], [1]),
    Spec.HoppingOffdiagonal( t2, "A", "A", [0], [2]),
  ]
  func1 = Generator.generatefast(uc, hopspec)

  hopembed = Embed.Hopping[Embed.embed(uc, hop) for hop in hopspec]
  func2 = Generator.generatefast(uc, hopembed)

  for kx in linspace(0.0, 2.0*pi, 16+1)
    out1 = zeros(Complex128, (1,1))
    out2 = zeros(Complex128, (1,1))
    func1([kx], out1)
    func2([kx], out2)
    ek = (t0
          + t1 * exp(1im * kx) + conj(t1) * exp(-1im * kx)
          + t2 * exp(2im * kx) + conj(t2) * exp(-2im * kx) )
    #@show out[1,1]
    #@show ek
    @test isapprox(out1[1,1], ek)
    @test isapprox(out2[1,1], ek)
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
  uc = UnitCell([b1 b3])
  addorbital!(uc, "A", carte2fract(uc, [0.0, 0.0]))
  addorbital!(uc, "B", carte2fract(uc, [0.0,-1.0]))
  #@show uc

  hopspecs = Spec.Hopping[]
  t1 = 1.0 + 0.1im
  for r in [a1, a2, a3]
    push!( hopspecs, Spec.hoppingbycarte(uc, t1, "A", "B", [0.0, 0.0], r) )
  end

  func1 = Generator.generatefast(uc, hopspecs)

  hopembed = Embed.Hopping[Embed.embed(uc, hop) for hop in hopspecs]
  func2 = Generator.generatefast(uc, hopembed)

  for kx in linspace(4.0, 4.0, 8)
    for ky in linspace(-4.0, 4.0, 8)
      k = [kx, ky]
      mat = begin
        phase = (exp(1im * dot(k, a1)) + exp(1im * dot(k, a2)) + exp(1im * dot(k, a3)) )
        T1 = t1 * phase
        T1c = conj(T1)
        [  0  T1 ;  T1c  0 ]
      end

      out1 = zeros(Complex128, (2, 2))
      out2 = zeros(Complex128, (2, 2))
      func1(k, out1)
      func2(k, out2)

      @test isapprox(out1, mat)
      @test isapprox(out2, mat)
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
  hopembed = Embed.Hopping[Embed.embedhopping(uc, hop) for hop in hopspec]

  func = Generator.generate(uc, hopembed)
  out = zeros(Complex128, (2,2))
  @show out
  func([0.0, 0.0], out)
  @show out

  out = zeros(Complex128, (2,2))
  @show out
  func([0.0001, 0.0002], out)
  @show out

  out = zeros(Complex128, (2,2))
  @show out
  func([0.0001, 0.0002], out)
  @show out

end

=#
