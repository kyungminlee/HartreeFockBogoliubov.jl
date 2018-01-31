module LinearizedGap

using HartreeFockBogoliubov


using ..Lattice
import ..Spec
import ..HFB

export linearizedpairingkernel

""" linearizedpairingkernel

Compute the kernel Γ of the linearized gap equation in the pairing channel
which is written as Δ = Γ⋅Δ
"""
function linearizedpairingkernel(
            solver::HFB.HFBSolver{O};
            tolerance::Float64 = sqrt(eps(Float64))
        ) where {O}

    unitcell = solver.hamiltonian.unitcell
    size = solver.size
    hfbcomputer = solver.hfbcomputer

    num_orbital = numorbital(unitcell)
    num_eigen = num_orbital
    num_pairing = length(hfbcomputer.Δ_registry)

    momentumgrid = Lattice.momentumgrid(unitcell, size)
    #@show momentumgrid
    num_momentum = length(momentumgrid)

    hoppinghamiltonian = HFB.makehoppingmatrix(hfbcomputer)

    # 1. Compute eigenenergies
    # Eigenstates of non-interacting Hamiltonian H0
    eigenenergies_plusk = zeros(Float64, (num_orbital, num_momentum))
    eigenstates_plusk = zeros(Complex128, (num_orbital, num_orbital, num_momentum))
    eigenenergies_minusk = zeros(Float64, (num_orbital, num_momentum))
    eigenstates_minusk = zeros(Complex128, (num_orbital, num_orbital, num_momentum))

    for (idx_momentum, momentum) in enumerate(momentumgrid)
        hk = hoppinghamiltonian(momentum)
        (w, v) = eig(Hermitian(hk))
        eigenenergies_plusk[:, idx_momentum] = w
        eigenstates_plusk[:, :, idx_momentum] = v

        hk = hoppinghamiltonian(-momentum)
        (w, v) = eig(Hermitian(hk))
        eigenenergies_minusk[:, idx_momentum] = w
        eigenstates_minusk[:, :, idx_momentum] = v
    end
    fermis_plusk = hfbcomputer.fermi.(eigenenergies_plusk)
    fermis_minusk = hfbcomputer.fermi.(eigenenergies_minusk)

    lindhards = zeros(Float64, (num_eigen, num_eigen, num_momentum))
    for (ik, momentum) in enumerate(momentumgrid)
        for n in 1:num_eigen, m in 1:num_eigen
            E1 = eigenenergies_plusk[n, ik]
            E2 = eigenenergies_minusk[m, ik]
            f1 = fermis_plusk[n, ik]
            f2 = fermis_minusk[m, ik]
            if abs(E1 + E2) > tolerance
                #@show ( 1.0 - f1 - f2 ) / ( E1 + E2 )
                lindhards[n, m, ik] = ( 1.0 - f1 - f2 ) / ( E1 + E2 )
            else
                T = hfbcomputer.temperature
                if abs(T) > tolerance
                    lindhards[n, m, ik] = 0.5 / (T * (1 + cosh(E1 / T)))
                else
                    warn("Lindhard: Temperature and energy difference both zero")
                    lindhards[n, m, ik] = 0.5 / (tolerance * (1 + cosh(E1 / tolerance)))
                end
            end
        end
    end

    ##@show hfbcomputer.t_registry
    kernel_pp = zeros(Complex128, (num_pairing, num_pairing))
    for (ibond3, (_, i3, j3, rij3, _)) in enumerate(hfbcomputer.Δ_registry)
        for (ibond1, (_, i1, j1, rij1, srcs)) in enumerate(hfbcomputer.Δ_registry)
            for (srcidx, v, neg) in srcs
                (_, i2, j2, rij2) = hfbcomputer.t_registry[srcidx]

                kernel_element = zero(Complex128)
                for (ik, momentum) in enumerate(momentumgrid)
                    phase1 = exp(1im * dot(momentum, -rij2 + rij3) )
                    phase2 = exp(1im * dot(momentum, -rij2 - rij3) )
                    for n in 1:num_eigen, m in 1:num_eigen
                        ψ1 = eigenstates_plusk[:, n, ik]
                        ψ2 = eigenstates_minusk[:, m, ik]

                        lindhard = lindhards[n, m, ik]

                        kernel_element += lindhard * (
                        ψ1[i2] * ψ2[j2] * conj(ψ1[i3] * ψ2[j3]) * phase1
                        - ψ1[i2] * ψ2[j2] * conj(ψ1[j3] * ψ2[i3]) * phase2 )
                    end
                end
                kernel_element /= length(momentumgrid)
                kernel_pp[ibond1, ibond3] += v * (neg ? -1 : 1) * kernel_element
            end
        end
    end
    return kernel_pp
end


end # module LinearizedGap
