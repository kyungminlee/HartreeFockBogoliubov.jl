
using ..Generator
using LinearAlgebra

"""
    chernnumber

Compute chern number of the band structure defined by the hoppings and the selected bands.

"""
function chernnumber(uc::UnitCell{O},
                     hoppings::AbstractVector{Hopping},
                     n1 ::Integer,
                     n2 ::Integer,
                     selectband::AbstractVector{<:Integer};
                     tol::Real=sqrt(eps(Float64))) where {O}
    squareuc = squarify(uc)
    hkgen = Generator.generatefast(squareuc, hoppings)
    norb = numorbital(squareuc)

    kgrid = momentumgrid(squareuc, [n1, n2])
    nk = length(kgrid)
    nsel = length(selectband)

    eigenvaluegrid = zeros(Float64, (nk, nsel))
    eigenvectorgrid = zeros(Complex{Float64}, (nk, norb, nsel))

    hk = zeros(Complex{Float64}, (norb, norb))
    for (i, k) in enumerate(kgrid)
        hk[:] = 0.0
        hkgen(k, hk)
        u, v = eig(Hermitian(0.5 * (hk + hk')))
        eigenvaluegrid[i, :]     = u[selectband]
        eigenvectorgrid[i, :, :] = v[:, selectband]
    end

    eigenvaluegrid  = reshape(eigenvaluegrid,  (n1, n2, nsel))
    eigenvectorgrid = reshape(eigenvectorgrid, (n1, n2, norb, nsel))

    flux = 0.0
    for i2 in 1:n2, i1 in 1:n1
        i1p = mod1(i1+1, n1)
        i2p = mod1(i2+1, n2)
        ψ1 = eigenvectorgrid[i1 , i2 , :, :]
        ψ2 = eigenvectorgrid[i1p, i2 , :, :]
        ψ3 = eigenvectorgrid[i1p, i2p, :, :]
        ψ4 = eigenvectorgrid[i1 , i2p, :, :]
        flux += angle(det(ψ1' * ψ2) * det(ψ2' * ψ3) * det(ψ3' * ψ4) * det(ψ4' * ψ1))
    end
    chn = flux / (2π)
    chnint = round(chn)
    if ! isapprox(chn, chnint; atol=tol)
        @warn( "Chern number should be almost integer")
    end
    return chn
end
