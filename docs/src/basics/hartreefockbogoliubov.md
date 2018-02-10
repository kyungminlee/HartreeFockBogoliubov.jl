# Hartree-Fock-Bogoliubov Theory

This document contains basics of Hartree-Fock-Bogoliubov Theory. We will follow [Goodman][Goodman80].

```math
  H =  \sum_{ij} T_{ij} c_{i}^{*} c_{j}
     + \frac{1}{4} \sum_{ijkl} V_{ijkl} c_{i}^{*} c_{j}^{*} c_{l} c_{k}
```
The finite-temperature properties of the system described by the above Hamiltonian can be computed using the density matrix
```math
\begin{align}
D &= \frac{Z} e^{-\beta H} \\
Z = \mathrm{Tr} e^{-\beta H}
\end{align}
```
In Hartree-Fock-Bogoliubov theory, we approximate the Hamiltonian ``H`` by a non-interacting Hamiltonian of "quasiparticles"
```math
H_{\text{HFB}} = E_0 + \sum_{n} E_{n} a_{n}^{*} a_{n}
```
The quasiparticle operator ``a`` is related to ``c`` in the following way:
```math
c_{i} = \sum_{n} \left( U_{in} a_{n} + V_{in} a_{n}^{*} \right)
```
The HFB density matrix is given by
```math
\begin{align}
  D_{\text{HFB}} &= \frac{1}{Z_{\text{HFB}}} e^{-\beta \sum_{n} E_{n} \hat{n}_{n} } \\
  Z_{\text{HFB}} &= \mathrm{Tr} e^{-\beta \sum_{n} E_{n} \hat{n}_{n} }
\end{align}
```
where
```math
\hat{n}_{n} = a_{n}^{*} a_{n}
```
is the quasiparticle number operator.

The HFB partition function can be calculated:
```math
Z_{\text{HFB}} = \prod_{n} \left( 1 + e^{-\beta E_{n}} \right)
```
and the density matrix is given in terms of the quasiparticle number operator
```math
D_{\text{HFB}} = Z_{\text{HFB}}^{-1}
\prod_{n}
\left[
  e^{-\beta E_{n}} \hat{n}_{n} + (1 - \hat{n}_{n})
\right]
= \prod_{n} \left[ f_{n} \hat{n}_{n} + (1-f_{n}) (1 - \hat{n}_{n}) \right]
```
where
```math
f_{n} = \frac{1}{1 + e^{\beta E_{n}}}
```
is the Fermi-Dirac distribution function for the nth quasiparticle

Single-quasparticle density matrix ``\overline{\rho}``  and pairing tensor ``\overline{t}``
```math
\begin{align}
\overline{\rho}_{ij}
  &= \left\langle a_{j}^{*} a_{i} \right\rangle
   = \mathrm{Tr} \left( D a_{j}^{*} a_{i} \right) \\
\overline{t}_{ij}
  &= \left\langle a_{j} a_{i} \right\rangle
   = \mathrm{Tr} \left( D a_{j} a_{i} \right)
\end{align}
```
Within HFB,
```math
\begin{align}
  \overline{\rho}_{ij} = \delta_{ij} f_{i} \\
  \overline{t}_{ij} = 0
\end{align}
```
The single-particle density matrix and pairing tensor are
```math
\begin{align}
\rho_{ij}
  &= \left\langle c_{j}^{*} c_{i} \right\rangle
  = \mathrm{Tr} \left( D c_{j}^{*} c_{i} \right) \\
t_{ij}
  &= \left\langle c_{j} c_{i} \right\rangle
  = \mathrm{Tr} \left( D c_{j} c_{i} \right)
\end{align}
```

```math
\begin{align}
  \rho &= U f U^{\dagger} + V (1-f) V^{\dagger} \\
  t    &= U f V^{\intercal} + V (1-f) U^{\intercal}
\end{align}
```
where ``f_{ij} = \delta_{ij} f_{i}``.



### Expectation Values
Wick's theorem
```math
\left\langle c_{i}^{*} c_{j}^{*} c_{l} c_{k} \right\rangle
=
\left\langle c_{i}^{*} c_{k} \right\rangle
\left\langle c_{j}^{*} c_{l} \right\rangle
- \left\langle c_{i}^{*} c_{l} \right\rangle
\left\langle c_{j}^{*} c_{k} \right\rangle
+ \left\langle c_{i}^{*} c_{j}^{*} \right\rangle
\left\langle c_{l} c_{k} \right\rangle
```

```math
\begin{align}
E &= \mathrm{tr}
\left[
  \left( T + \frac{1}{2} \Gamma \right) \rho + \frac{1}{2} \Delta t^{\dagger}
\right] \\
S &= - k_B \sum_{i} \left[ f_{i} \ln f_{i} + (1-f_{i}) \ln (1-f_i)\right] \\
N &= \mathrm{tr} \rho
\end{align}
```

```math
\begin{align}
\Gamma_{ij} = \sum_{kl} V_{ikjl} \rho_{lk} \\
\Delta_{ij} = \frac{1}{2} \sum_{kl} V_{ijkl} t_{kl} \\
\end{align}
```

The grand potential
```math
\Omega
= \sum_{ij} (T - \mu)_{ij} \rho_{ji}
+ \frac{1}{2} \sum_{ijkl} V_{ijkl} \rho_{lj} \rho_{ki}
+ \frac{1}{4} \sum_{ijkl} V_{ijkl} t_{ij}^{*} t_{kl}
+ k_B T \sum_{i} \left[ f_{i} \ln f_{i} + (1-f_{i}) \ln (1-f_i)\right]
```


----





[Goodman80]: http://escholarship.org/uc/item/3xm630cr "Alan Goodman, Finite-Temperature Hartree-Fock-Bogoliubov Theory, LBNL Paper LBL-11151 (1980)"
