# Momentum Space Formulation

Hopping Elements

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


Expectation Values

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


