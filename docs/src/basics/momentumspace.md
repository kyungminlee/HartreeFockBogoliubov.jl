# Momentum Space Formulation

Hopping Elements

```math
\sum_{\bfR}
    c_{\alpha}^{*}(\bfR + \boldsymbol{\rho}_\alpha)
    c_{\beta}(\bfR + \boldsymbol{\rho}_\beta)
  =
    \sum_{\bfk}
    e^{i \bfk \cdot (\boldsymbol{\rho}_\beta - \boldsymbol{\rho}_\alpha)}
    c_{\alpha}^{*}(\bfk)
    c_{\beta}(\bfk)
```

Pairing Elements
```math
\sum_{\bfR}
    c_{\alpha}^{*}(\bfR + \boldsymbol{\rho}_\alpha)
    c_{\beta}^{*}(\bfR + \boldsymbol{\rho}_\beta)
  =
    \sum_{\bfk}
    e^{i \bfk \cdot (\boldsymbol{\rho}_\beta - \boldsymbol{\rho}_\alpha)}
    c_{\alpha}^{*}(\bfk)
    c_{\beta}^{*}(-\bfk)
```


```math
  \sum_{\bfR}
    c_{\alpha}(\bfR + \boldsymbol{\rho}_\alpha)
    c_{\beta}(\bfR + \boldsymbol{\rho}_\beta) 
  =
    \sum_{\bfk}
    e^{i \bfk \cdot (\boldsymbol{\rho}_\beta - \boldsymbol{\rho}_\alpha)}
    c_{\alpha}(-\bfk)
    c_{\beta}(\bfk)
```


Expectation Values

```math
  \rho_{\alpha \beta}
  =
      \frac{1}{N}
      \sum_{\bfk}
      e^{-i \bfk \cdot ( \boldsymbol{\rho}_{\beta} - \boldsymbol{\rho}_{\alpha} ) }
      \sum_{n}
       f(\epsilon_n) U_{\alpha n} U_{\beta n}^{*}
```


```math
  t_{\alpha\beta}
  =
      \frac{1}{N}
      \sum_{\bfk_1, \bfk_2}
      e^{-i \bfk \cdot ( \boldsymbol{\rho}_{\beta}
                        -\boldsymbol{\rho}_{\alpha} ) }
      \sum_{n}
        f(\epsilon_n )
        U_{\alpha n}
        V_{\beta n}^{*}
```


