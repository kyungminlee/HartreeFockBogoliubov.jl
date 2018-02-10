# Hartree-Fock-Bogoliubov Theory

This document contains basics of Hartree-Fock-Bogoliubov Theory. We will follow [Goodman][Goodman80].

```math
  H =  \sum_{ij} T_{ij} c_{i}^{\dagger} c_{j}
     + \frac{1}{4} \sum_{ijkl} V_{ijkl} c_{i}^{\dagger} c_{j}^{\dagger} c_{l} c_{k}
```

```math
H_{\text{HFB}} = E_0 + \sum_{i} E_{n} a_{n}^{\dagger} a_{n}
```

```math
c_{i} = \sum_{n} \left( U_{in} a_{n} + V_{in} a_{n}^{*} \right)
```
TODO: is this correct?


```math
\begin{align}
  D_{\text{HFB}} &= \frac{1}{Z_{\text{HFB}}} e^{-\beta \sum_{i} E_{i} \hat{n}_{i} } \\
  Z_{\text{HFB}} &= \mathrm{Tr} e^{-\beta \sum_{i} E_{i} \hat{n}_{i} }
\end{align}
```

```math
\hat{n}_{i} = a_{i}^{\dagger} a_{i}
```


```math
Z_{\text{HFB}} = \prod_{i} \left( 1 + e^{-\beta E_{i}} \right)
```


```math
D_{\text{HFB}} = Z_{\text{HFB}}^{-1} \prod_{i}
\left[
  e^{-\beta E_{i}} \hat{n}_{i} + (1 - \hat{n}_i)
\right]
= \prod_{i} \left[ f_{i} \hat{n}_{i} + (1-f_{i}) (1 - \hat{n}_{i}) \right]
```
where
```math
f_{i} = \frac{1}{1 + e^{\beta E_{i}}}
```


Single-quasparticle density matrix ``\overline{\rho}``  and pairing tensor ``\overline{t}``
```math
\begin{align}
\overline{\rho}_{ij}
  &= \left\langle a_{j}^{\dagger} a_{i} \right\rangle
   = \mathrm{Tr} \left( D a_{j}^{\dagger} a_{i} \right) \\
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
  &= \left\langle c_{j}^{\dagger} c_{i} \right\rangle
  = \mathrm{Tr} \left( D c_{j}^{\dagger} c_{i} \right) \\
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
\left\langle c_{i}^{\dagger} c_{j}^{\dagger} c_{l} c_{k} \right\rangle
=
\left\langle c_{i}^{\dagger} c_{k} \right\rangle
\left\langle c_{j}^{\dagger} c_{l} \right\rangle
- \left\langle c_{i}^{\dagger} c_{l} \right\rangle
\left\langle c_{j}^{\dagger} c_{k} \right\rangle
+ \left\langle c_{i}^{\dagger} c_{j}^{\dagger} \right\rangle
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
\Omega = \sum_{ij} (T - \mu)_{ij} \rho_{ji} + \frac{1}{2} \sum_{ijkl} V_{ijkl} \rho_{l j} \rho_{ki}
+ \frac{1}{4} \sum_{ijkl} V_{ijkl} t_{ij}^{*} t_{kl}
+ k_B T \sum_{i} \left[ f_{i} \ln f_{i} + (1-f_{i}) \ln (1-f_i)\right]
```
----








[Goodman80]: http://escholarship.org/uc/item/3xm630cr "Alan Goodman, Finite-Temperature Hartree-Fock-Bogoliubov Theory, LBNL Paper LBL-11151 (1980)"
