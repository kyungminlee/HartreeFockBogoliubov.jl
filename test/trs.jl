using LinearAlgebra
using DataStructures
using HartreeFockBogoliubov.Topology

function chopsmall(arr ::AbstractArray{C}; tolerance=sqrt(eps(Float64))) where {C<:Real}
    out = similar(arr)
    for (i, x) in enumerate(arr)
        out[i] = abs(x) > tolerance ? x : 0
    end
    return out
end

function chopsmall(arr ::AbstractArray{C}; tolerance=sqrt(eps(Float64))) where {C<:Complex}
    out = similar(arr)
    for (i, x) in enumerate(arr)
        rv = real(x)
        rv = abs(rv) > tolerance ? rv : 0
        iv = imag(x)
        iv = abs(iv) > tolerance ? iv : 0
        out[i] = rv + im * iv
    end
    return out
end

function main()
    srand(1)

    uT = [  0  0  1  0  0  0  0  0  ;
            0  0  0  1  0  0  0  0  ;
           -1  0  0  0  0  0  0  0  ;
            0 -1  0  0  0  0  0  0  ;
            0  0  0  0  0  0 -1  0  ;
            0  0  0  0  0  0  0 -1  ;
            0  0  0  0  1  0  0  0  ;
            0  0  0  0  0  1  0  0  ]

    @show isvalidtimereversalmatrix(uT)

    v1 = zeros(ComplexF64, 8, 8)

    # TODO ALSO CHECK  COMPLEX
    for idxpair in 1:4
        vec = rand(ComplexF64, 8)
        for idxprev in 1:(2*idxpair-1)
            vec -= v1[:, idxprev] * dot(v1[:, idxprev], vec)
        end
        normalize!(vec)
        v1[:, 2*idxpair-1] = vec
        v1[:, 2idxpair] = uT * conj(vec)
    end

    w1 = [1, 1, 2, 2, 3, 3, 3, 3]
    h = v1 * Diagonal(w1) * v1'
    h = 0.5 * (h + h')

    println("h")
    display(chopsmall(h))
    println()


    @show istimereversal(h, uT)
    (w2, v2) = eig(Hermitian(h))

    println("v2")
    display( chopsmall(v2) )
    println()

    v3 = copy(v2)
    kramerpairup!(w2, v3, uT)


    @show w1
    @show w2

    println("v1")
    display( chopsmall(v1) )
    println()

    println("v2")
    display( chopsmall(v2) )
    println()

    println("v3")
    display( chopsmall(v3) )
    println()

    println("v1 * v1'")
    display( chopsmall(v1 * v1') )
    println()
    println("v1' * v1")
    display( chopsmall(v1' * v1) )
    println()

    println("v2 * v2'")
    display( chopsmall(v2 * v2') )
    println()
    println("v2' * v2")
    display( chopsmall(v2' * v2) )
    println()

    println("v3 * v3'")
    display( chopsmall(v3 * v3') )
    println()
    println("v3' * v3")
    display( chopsmall(v3' * v3) )
    println()

    println("v2' * v3")
    display(chopsmall(v2' * v3))
    println()

    println("CANONICAL?")

    println("transpose(v1) * uT * v1")
    display(chopsmall(transpose(v1) * uT * v1))
    println()

    println("transpose(v2) * uT * v2")
    display(chopsmall(transpose(v2) * uT * v2))
    println()

    println("transpose(v3) * uT * v3")
    display(chopsmall(transpose(v3) * uT * v3))
    println()

end


function main2()
  srand(1)

  uT = [ 0  1  0  0 ;
        -1  0  0  0 ;
         0  0  0 -1 ;
         0  0  1  0 ]

  v = zeros(Float64, 4, 4)

  let
    vec = rand(Float64, 4)
    normalize!(vec)
    v[:, 1] = vec
    v[:, 2] = uT * conj(vec)
  end

  let
    vec = rand(Float64, 4)
    vec -= v[:, 1] * dot(v[:,1], vec) + v[:, 2] * dot(v[:, 2], vec)
    normalize!(vec)
    v[:, 3] = vec
    v[:, 4] = uT * conj(vec)
  end

  h = v * v' * 0.2

  @show h

  @assert( istimereversal(h, uT) )

  (w, v) = eig(Hermitian(h))

  @show w
  @show "BEFORE:", v
  @show chopsmall( v * v' )
  @show chopsmall( v' * v )

  res = kramerpairupblock!(v, uT)

  @show "AFTER:", v
  @show chopsmall( v * v' )
  @show chopsmall( v' * v )
  #@show res
  #@show res * res'
  #@show res' * res

  #=
  h = rand(Float64, 4, 4)
  h = 0.5 * (h + h')

  # want h - h2 = dh = 0

  for run in 1:100
    h2 = conj(u * h * u')
    dh = h - h2
    ndh = norm(dh)
    h -= 0.5 * dh
    if isapprox(ndh, 0)
      break
    end
  end

  for run in 1:4
    h2 = conj(u * h * u')
    dh = h - h2
    ndh = norm(dh)
    h -= 0.5 * dh
    @show ndh
    if isapprox(ndh, 0)
      break
    end
  end

  h = 0.5 * (h + h')
  @show norm( h - conj(u * h * u') )

  @show h

  # Let's hope that h is time-reversal.

  (w, v) = eig(Hermitian(h))
  @show w
  @show v


  @show kramer(u, v)
  =#

  nothing
end

main()
