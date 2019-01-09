# Momentum Space Formulation

We can make use of the translation symmetry of the system by working in momentum space.

## Self-consistency Loop

### Hamiltonian in Momentum Space

Translation in the hopping/pairing Hamiltonian implies that the hopping/pairing terms are summed over all the (Bravais) lattice vectors. Upon Fourier transform, the hopping terms become
```math
\sum_{\mathbf{R}}
    c_{\alpha}^{*}(\mathbf{R} + \boldsymbol{\rho}_\alpha)
    c_{\beta}(\mathbf{R} + \boldsymbol{\rho}_\beta)
  =
    \sum_{\mathbf{k}}
    e^{i \mathbf{k} \cdot (\boldsymbol{\rho}_\beta - \boldsymbol{\rho}_\alpha)}
    c_{\alpha}^{*}(\mathbf{k})
    c_{\beta}(\mathbf{k})
```
and the pairing terms become
```math
\sum_{\mathbf{R}}
    c_{\alpha}^{*}(\mathbf{R} + \boldsymbol{\rho}_\alpha)
    c_{\beta}^{*}(\mathbf{R} + \boldsymbol{\rho}_\beta)
  =
    \sum_{\mathbf{k}}
    e^{i \mathbf{k} \cdot (\boldsymbol{\rho}_\beta - \boldsymbol{\rho}_\alpha)}
    c_{\alpha}^{*}(\mathbf{k})
    c_{\beta}^{*}(-\mathbf{k})
```
and
```math
\sum_{\mathbf{R}}
    c_{\alpha}(\mathbf{R} + \boldsymbol{\rho}_\alpha)
    c_{\beta}(\mathbf{R} + \boldsymbol{\rho}_\beta)
  =
    \sum_{\mathbf{k}}
    e^{i \mathbf{k} \cdot (\boldsymbol{\rho}_\beta - \boldsymbol{\rho}_\alpha)}
    c_{\alpha}(-\mathbf{k})
    c_{\beta}(\mathbf{k})
```


## Expectation Values in Momentum Space

If the Hamiltonian has translation symmetry, the expectation value of the hopping and pairing amplitudes (which are assumed to be the same on every unitcell, of course) can be written as a sum of the expectation values in each momentum sector. The hopping amplitude (or density if diagonal) is
```math
\rho_{\alpha \beta}
  =
      \frac{1}{N}
      \sum_{\mathbf{k}}
      e^{-i \mathbf{k} \cdot ( \boldsymbol{\rho}_{\beta} - \boldsymbol{\rho}_{\alpha} ) }
      \sum_{n}
       f(\epsilon_n) U_{\alpha n} U_{\beta n}^{*}
```
whereas the pairing amplitude is
```math
t_{\alpha\beta}
  =
      \frac{1}{N}
      \sum_{\mathbf{k}_\alpha, \mathbf{k}_\beta}
      e^{-i \mathbf{k} \cdot ( \boldsymbol{\rho}_{\beta}
                        -\boldsymbol{\rho}_{\alpha} ) }
      \sum_{n}
        f(\epsilon_n )
        U_{\alpha n}
        V_{\beta n}^{*}
```


## Free Energy

The Hartree-Fock-Bogoliubov can be solved using only the self-consistent loop. There is, however, always a danger of being trapped in a local optimum. In order to get around this issue, we can compute the Hartree-Fock-Bogoliubov free energy of the solution, and compare them to make sure optimal solution has been found.

### Grand Potential

The energy of a HFB solution can be written as
```math
E = \mathrm{tr} \left[ \left( T + \frac{1}{2} \Gamma \right) \rho + \frac{1}{2} \Delta t^{\dagger} \right]
  = E_{T} + E_{\Gamma} + E_{\Delta}
```
where
```math
\begin{align}
E_{T} &= \mathrm{tr} \left( T \rho \right) \\
E_{\Gamma} &= \frac{1}{2} \mathrm{tr} \left( \Gamma \rho        \right) \\
E_{\Delta} &= \frac{1}{2} \mathrm{tr} \left( \Delta t^{\dagger} \right)
\end{align}
```

In momentum space
```math
\begin{align}
T_{ij}^{\mathbf{k}} &= T_{ij} e^{i \mathbf{k} \cdot \boldsymbol{\rho}_{ij}} \\
\rho_{ij} &= \frac{1}{N} \sum_{\mathbf{k}} \rho_{ij}^{\mathbf{k}} e^{-i \mathbf{k} \cdot \boldsymbol{\rho}_{ij}}
\end{align}
```
Thus the expectation value of the kinetic term is
```math
\begin{align}
E_{T}
  &= \sum_{ij} T_{ij} \rho_{ji}
   = \sum_{ij} T_{ij}
     \frac{1}{N} \sum_{\mathbf{k}} \rho_{ji}^{\mathbf{k}} e^{-i \mathbf{k} \cdot \boldsymbol{\rho}_{ji} }
   = \frac{1}{N} \sum_{\mathbf{k}} \sum_{ij} T_{ij} e^{i \mathbf{k} \cdot \boldsymbol{\rho}_{ij} }
      \rho_{ji}^{\mathbf{k}}
   = \frac{1}{N} \sum_{ij}
      \sum_{\mathbf{k}}
      T_{ij}^{\mathbf{k}}
      \rho_{ji}^{\mathbf{k}}
\end{align}
```
Similarly for Γ and Δ:
```math
\begin{align}
E_{\Gamma}
  &= \frac{1}{2} \frac{1}{N} \sum_{ij}
    \sum_{\mathbf{k}}
    \Gamma_{ij}^{\mathbf{k}}
    \rho_{ji}^{\mathbf{k}} \\
E_{\Delta}
  &= \frac{1}{2} \frac{1}{N} \sum_{ij}
    \sum_{\mathbf{k}}
    \Delta_{ij}^{\mathbf{k}}
    {t_{ij}^{\mathbf{k}} }^{*} \\
\end{align}
```
One thing to be careful is that, when using these formula, (Γ, Δ) need to be consistent with (ρ, t); in other words, one needs to use the Γ and Δ calculated from the
ρ and t. If, on the other hand, one uses (ρ, t) calculated from the HFB Hamiltonian generated with (Γ, Δ), the calculated grand potential may, in fact, be lower than that of the true (HFB) ground state. Of course, this is not an issue once self-consistency is reached.
