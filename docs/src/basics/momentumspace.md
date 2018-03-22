# Momentum Space Formulation

### Hopping Elements

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

Pairing Elements
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


### Expectation Values

```math
\rho_{\alpha \beta}
  =
      \frac{1}{N}
      \sum_{\mathbf{k}}
      e^{-i \mathbf{k} \cdot ( \boldsymbol{\rho}_{\beta} - \boldsymbol{\rho}_{\alpha} ) }
      \sum_{n}
       f(\epsilon_n) U_{\alpha n} U_{\beta n}^{*}
```


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



### Grand Potential

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
T_{ij}^{\mathbf{k}} &= T_{ij} e^{i \mathbf{k} \cdot \boldsymbol{\rho}_{ij}}
\rho_{ij} = \frac{1}{N} \sum_{\mathbf{k}} \rho_{ij}^{\mathbf{k}} e^{-i \mathbf{k} \cdot \boldsymbol{\rho}_{ij}}
\end{align}
```

Thus
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
Here however, Γ and Δ need to be consistent with ρ and t.
