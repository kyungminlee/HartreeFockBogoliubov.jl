
# Hartree-Fock-Bogoliubov Decomposition

We can write the Hamiltonian in the way that makes the Hermiticity manifest, taking into account duplicates properly
``` math
\begin{align}
H &= \sum_{ij} T_{ij} c_{i}^{*} c_{j}
 + \frac{1}{4} \sum_{ijkl} c_{i}^{*} c_{j}^{*} c_{l} c_{k}
\\
  &= \sum_{i} T_{ii} c_{i}^{*} c_{i}
    + \sum_{i \lt j} (T_{ij} c_{i}^{*} c_{j} + T_{ij}^{*} c_{j}^{*} c_{i}) \\
  &\quad
  + \frac{1}{4} \sum_{i < j}
  \left[
      V_{ijij} c_{i}^{*} c_{j}^{*} c_{j} c_{i}
    - V_{ijij} c_{j}^{*} c_{i}^{*} c_{j} c_{i}
    - V_{ijij} c_{i}^{*} c_{j}^{*} c_{i} c_{j}
    + V_{ijij} c_{j}^{*} c_{i}^{*} c_{i} c_{j}
  \right]
  \\
  &\quad
  + \frac{1}{4} \sum_{\substack{(i<j), (k<l) \\ \text{ and } \\ (i,j) < (k,l) }}
  \left[
      V_{ijkl}     c_{i}^{*} c_{j}^{*} c_{l} c_{k}
    - V_{ijkl}     c_{j}^{*} c_{i}^{*} c_{l} c_{k}
    - V_{ijkl}     c_{i}^{*} c_{j}^{*} c_{k} c_{l}
    + V_{ijkl}     c_{j}^{*} c_{i}^{*} c_{k} c_{l} \right.\\
  &\qquad \left.
    + V_{ijkl}^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}
    - V_{ijkl}^{*} c_{l}^{*} c_{k}^{*} c_{j} c_{i}
    - V_{ijkl}^{*} c_{k}^{*} c_{l}^{*} c_{i} c_{j}
    + V_{ijkl}^{*} c_{l}^{*} c_{k}^{*} c_{i} c_{j}
  \right]
\\
&= \sum_{i} T_{ii} c_{i}^{*} c_{i}
  + \sum_{i \lt j} (T_{ij} c_{i}^{*} c_{j} + T_{ij}^{*} c_{j}^{*} c_{i}) \\
&\quad
  + \sum_{i < j}
  \left[
      V_{ijij} c_{i}^{*} c_{j}^{*} c_{j} c_{i}
  \right]
  + \sum_{\substack{(i<j), (k<l) \\ \text{ and } \\ (i,j) < (k,l) }}
  \left[
      V_{ijkl}     c_{i}^{*} c_{j}^{*} c_{l} c_{k}
    + V_{ijkl}^{*} c_{k}^{*} c_{l}^{*} c_{j} c_{i}
  \right]
\end{align}
```

```math
\begin{align}
  \Gamma_{ij} &= \sum_{kl} V_{ikjl} \rho_{lk} \\
  \Gamma_{ik} &= \sum_{jl} V_{ijkl} \rho_{lj} \\
  \Delta_{ij} &= \frac{1}{2} \sum_{kl} V_{ijkl} t_{kl}
\end{align}
```


## Diagonal Interaction

Let us first consider the *diagonal* interaction of the form
```math
  V c_{1}^{*} c_{2}^{*} c_{2} c_{1} ,
```
which is the Hermitian conjugate of itself. Following [Goodman][Goodman80], this interaction decomposes into the following mean fields in the particle-hole and particle-particle channels:

| Interaction  | PH Channel       | PP Channel         |
|:------------ |:---------------- |:------------------ |
| V(1212) =  V | Γ(11) =  V ρ(22) | Δ(12) =  ½ V t(12) |
| V(1221) = -V | Γ(12) = -V ρ(12) | Δ(12) = -½ V t(21) |
| V(2112) = -V | Γ(21) = -V ρ(21) | Δ(21) = -½ V t(12) |
| V(2121) =  V | Γ(22) =  V ρ(11) | Δ(21) =  ½ V t(21) |

Here, the equal sign '=' does not indicate equality; it refers to the contribution of the right hand side to the mean field on the left hand side. Some of the mean fields are redundant, as required by the Hermiticity of the Hamiltonian. Overall, we end up with
```math
\begin{align}
\Gamma_{11} &=  V \rho_{22} \\
\Gamma_{12} &= -V \rho_{12} \\
\Gamma_{22} &=  V \rho_{11}
\end{align}
```
in the particle-hole channel, and
```math
\begin{align}
\Delta_{12} =  V t_{12}
\end{align}
```
in the particle-particle channel.


## Offdiagonal Interaction

Now let us consider the *offdiagonal* interaction term of the form
```math
  V c_{1}^{*} c_{2}^{*} c_{4} c_{3}.
```
Given this term in the interaction, it is implied that its Hermitian conjugate
```math
  V^{*} c_{3}^{*} c_{4}^{*} c_{2} c_{1}.
```
is also included in the Hamiltonian. The two terms (the explicit term and its Hermitian conjugate) decompose into

| Interaction   | PH Channel       | PP Channel         |
|:------------- |:---------------- |:------------------ |
| V(1234) =  V  | Γ(13) =  V ρ(42) | Δ(12) =  ½ V t(34) |
| V(1243) = -V  | Γ(14) = -V ρ(32) | Δ(12) = -½ V t(43) |
| V(2134) = -V  | Γ(23) = -V ρ(41) | Δ(21) = -½ V t(34) |
| V(2143) =  V  | Γ(24) =  V ρ(31) | Δ(21) =  ½ V t(43) |
|               |                  |                    |
| V(3412) =  V* | Γ(31) =  V ρ(24) | Δ(34) =  ½ V t(12) |
| V(3421) = -V* | Γ(32) = -V ρ(14) | Δ(34) = -½ V t(21) |
| V(4312) = -V* | Γ(41) = -V ρ(23) | Δ(43) = -½ V t(12) |
| V(4321) =  V* | Γ(42) =  V ρ(13) | Δ(43) =  ½ V t(21) |

As we have mentioned above for the diagonal interaction terms, this table contains redundant mean fields. Overall. we have
```math
\begin{align}
\Gamma_{13} &=  V  \rho_{24}^* \\
\Gamma_{14} &= -V  \rho_{23}^* \\
\Gamma_{23} &= -V  \rho_{14}^* \\
\Gamma_{24} &=  V  \rho_{13}^*
\end{align}
```
in the particle-hole channel and
```math
\begin{align}
\Delta_{12} &=  V   t_{34} \\
\Delta_{34} &=  V^* t_{12}
\end{align}
```
in the particle-particle channel.
