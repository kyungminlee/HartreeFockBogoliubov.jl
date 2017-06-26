# Hartree-Fock-Bogoliubov Theory

This document contains basics of Hartree-Fock-Bogoliubov Theory.

```math
  H =  \sum_{ij} T_{ij} c_{i}^{\dagger} c_{j}
     + \frac{1}{4} \sum_{ijkl} V_{ijkl} c_{i}^{\dagger} c_{j}^{\dagger} c_{l} c_{k}
```

```math
  \Gamma_{ij} = \sum_{kl} V_{ikjl} \rho_{lk}
```

```math
  \Gamma_{ik} = \sum_{kl} V_{ijkl} \rho_{lj}
```

```math
  \Delta_{ij} = \frac{1}{2} \sum_{kl} V_{ijkl} t_{kl}
```



## Hartree-Fock-Bogoliubov Decomposition


### Diagonal Interaction

```math
  V c_{1}^{\dagger} c_{2}^{\dagger} c_{2} c_{1}
```

| Interaction       | Γ                              | Δ                                       |
|:----------------- |:------------------------------ |:--------------------------------------- |
| ``V_{1212} =  V`` | ``\Gamma_{11} =  V \rho_{22}`` | ``\Delta_{12} =  \frac{1}{2} V t_{12}`` |
| ``V_{1221} = -V`` | ``\Gamma_{12} = -V \rho_{12}`` | ``\Delta_{12} = -\frac{1}{2} V t_{21}`` |
| ``V_{2112} = -V`` | ``\Gamma_{21} = -V \rho_{21}`` | ``\Delta_{21} = -\frac{1}{2} V t_{12}`` |
| ``V_{2121} =  V`` | ``\Gamma_{22} =  V \rho_{11}`` | ``\Delta_{21} =  \frac{1}{2} V t_{21}`` |

```math
\begin{aligned}
\Gamma_{11} =  V \rho_{22} \\
\Gamma_{12} = -V \rho_{12} \\
\Gamma_{22} =  V \rho_{11}
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
| ``V_{1234} =  V  `` | ``\Gamma{13} =  V \rho_{42}`` | ``\Delta_{12} =  \frac{1}{2} V t_{34}`` |
| ``V_{1243} = -V  `` | ``\Gamma{14} = -V \rho_{32}`` | ``\Delta_{12} = -\frac{1}{2} V t_{43}`` |
| ``V_{2134} = -V  `` | ``\Gamma{23} = -V \rho_{41}`` | ``\Delta_{21} = -\frac{1}{2} V t_{34}`` |
| ``V_{2143} =  V  `` | ``\Gamma{24} =  V \rho_{31}`` | ``\Delta_{21} =  \frac{1}{2} V t_{43}`` |
|:------------------- |:----------------------------- |:--------------------------------------- |
| ``V_{3412} =  V^*`` | ``\Gamma{31} =  V \rho_{24}`` | ``\Delta_{34} =  \frac{1}{2} V t_{12}`` |
| ``V_{3421} = -V^*`` | ``\Gamma{32} = -V \rho_{14}`` | ``\Delta_{34} = -\frac{1}{2} V t_{21}`` |
| ``V_{4312} = -V^*`` | ``\Gamma{41} = -V \rho_{23}`` | ``\Delta_{43} = -\frac{1}{2} V t_{12}`` |
| ``V_{4321} =  V^*`` | ``\Gamma{42} =  V \rho_{13}`` | ``\Delta_{43} =  \frac{1}{2} V t_{21}`` |


```math
\begin{aligned}
\Gamma_{13} =  V  \rho_{24}^* \\
\Gamma_{14} = -V  \rho_{23}^* \\
\Gamma_{23} = -V  \rho_{14}^* \\
\Gamma_{24} =  V  \rho_{13}^*
\end{aligned}
```
```math
\begin{aligned}
Δ_{12} =  V   t_{34} \\
Δ_{34} =  V^* t_{12}
\end{aligned}
```
