# Hartree-Fock-Bogoliubov Theory

This document contains basics of Hartree-Fock-Bogoliubov Theory.

```math
  H =  \sum_{ij} T_{ij} c_{i}^{\dagger} c_{j}
     + \frac{1}{4} \sum_{ijkl} V_{ijkl} c_{i}^{\dagger} c_{j}^{\dagger} c_{l} c_{k}
```

```math
\begin{aligned}
  \Gamma_{ij} &= \sum_{kl} V_{ikjl} \rho_{lk} \\
  \Gamma_{ik} &= \sum_{kl} V_{ijkl} \rho_{lj} \\
  \Delta_{ij} &= \frac{1}{2} \sum_{kl} V_{ijkl} t_{kl}
\end{aligned}
```



## Hartree-Fock-Bogoliubov Decomposition


### Diagonal Interaction

```math
  V c_{1}^{\dagger} c_{2}^{\dagger} c_{2} c_{1}
```

Inline Equation: ``V_{1212} =  V``

| Interaction       | Γ                              | Δ                                       |
|:----------------- |:------------------------------ |:--------------------------------------- |
| V(1212) =  V | Γ(11) =  V ρ(22) | Δ(12) =  ½ V t(12) |
| V(1221) = -V | Γ(12) = -V ρ(12) | Δ(12) = -½ V t(21) |
| V(2112) = -V | Γ(21) = -V ρ(21) | Δ(21) = -½ V t(12) |
| V(2121) =  V | Γ(22) =  V ρ(11) | Δ(21) =  ½ V t(21) |

```math
\begin{aligned}
\Gamma_{11} &=  V \rho_{22} \\
\Gamma_{12} &= -V \rho_{12} \\
\Gamma_{22} &=  V \rho_{11}
\end{aligned}
```

```math
\begin{aligned}
Δ_{12} =  V t_{12}
\end{aligned}
```

### Offdiagonal Interaction

```math
  V c_{1}^{\dagger} c_{2}^{\dagger} c_{4} c_{3}
```

| Interaction         | Γ                             | Δ                                       |
|:------------------- |:----------------------------- |:--------------------------------------- |
| V(1234) =  V   | Γ(13) =  V ρ(42) | Δ(12) =  ½ V t(34) |
| V(1243) = -V   | Γ(14) = -V ρ(32) | Δ(12) = -½ V t(43) |
| V(2134) = -V   | Γ(23) = -V ρ(41) | Δ(21) = -½ V t(34) |
| V(2143) =  V   | Γ(24) =  V ρ(31) | Δ(21) =  ½ V t(43) |
|:------------------- |:----------------------------- |:--------------------------------------- |
| V(3412) =  V* | Γ(31) =  V ρ(24) | Δ(34) =  ½ V t(12) |
| V(3421) = -V* | Γ(32) = -V ρ(14) | Δ(34) = -½ V t(21) |
| V(4312) = -V* | Γ(41) = -V ρ(23) | Δ(43) = -½ V t(12) |
| V(4321) =  V* | Γ(42) =  V ρ(13) | Δ(43) =  ½ V t(21) |


```math
\begin{aligned}
\Gamma_{13} &=  V  \rho_{24}^* \\
\Gamma_{14} &= -V  \rho_{23}^* \\
\Gamma_{23} &= -V  \rho_{14}^* \\
\Gamma_{24} &=  V  \rho_{13}^*
\end{aligned}
```
```math
\begin{aligned}
Δ_{12} =  V   t_{34} \\
Δ_{34} =  V^* t_{12}
\end{aligned}
```
