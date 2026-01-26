# FAIRR : Functional Approximation of Impulse Responses in R

`FAIRR` is an R package designed to estimate symmetric and asymmetric impulse response functions using the FAIR methodology (R. Barnichon, C. Matthes, 2018, Journal of Monetary Economics), using Gaussian basis functions :

$$\psi_h = \sum^K_{k=1} a_k \exp\left(\frac{(h-b)^2}{c^2}\right)$$

Or exponential basis functions :

$$\psi_h = a \exp(-\ell h)$$

The estimation is performed via Stan and users must therefore have a C++ compiler installed. If you encounter errors related to "Rtools" or "C++14", please ensure your toolchain is properly configured for Stan by following the [RStan Getting Started guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started). \## Update **1.1.0** : Added support for exponential basis functions. Changed the function names : Gaussian basis FAIR_G(K) can be estimated with `varfairG()` if you want to get orthogonal shocks from a VAR(p) estimated automatically and with `fairG()` if you already have a series of exogenous shocks. For exponential decay basis function FAIR_E(1), use `fairE()`, which requires shocks to already be exogenous or orthogonal.

## Update
**1.1.0** : Added support for exponential basis functions. Changed the function names : Gaussian basis FAIR_G(K) can be estimated with `varfairG()` if you want to get orthogonal shocks from a VAR(p) estimated automatically and with `fairG()` if you already have a series of exogenous shocks. For exponential decay basis function FAIR_E(1), use `fairE()`, which requires shocks to already be exogenous or orthogonal.

**1.1.1** : Fixed documentation wording to include the exponential basis function.


Plan for **1.2.0** (non-binding) : Adding up to K = 4 on Gaussian basis IRF's, adding a Gamma distribution basis function.

## Installation

```         
devtools::install_github("gauthiersr/FAIRR")
library(FAIRR)
```

## References

**Barnichon, R., & Matthes, C. (2018).** "Functional Approximation of Impulse Responses." *Journal of Monetary Economics*, 99, 41-55.

**Stan Development Team.** "RStan: the R interface to Stan." [mc-stan.org](https://mc-stan.org)
