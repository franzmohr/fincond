# fincond

Financial conditions indices from the time varying parameter FAVAR of
Koop and Korobilis (2014).

Model objects follow the conventions of
[bvartools](https://github.com/franzmohr/bvartools) and use its generics, so the
workflow is the one that package and [dfmtools](https://github.com/franzmohr/dfmtools)
use:

``` r
library("fincond")

data("koop")

model <- create_fcimodel(y = koop$macro[[96]], x = koop$financial, p = 4, n = 1)
model <- add_priors(model)
model <- add_posterior_coefficients(model)

index <- fci(model, sign_on = "vix")
plot(model)
```

`sign_on` fixes the orientation of the index, which is otherwise arbitrary: the
sign of a factor is not identified, so it has to be chosen rather than estimated.
Pinning it to the VIX makes a high value mean stressed conditions. On the paper's
own data the index peaks in 2008 Q4, ahead of 1982, 1974 and 1981.

`factor_loadings(model)` returns the exposure of each financial series to the
index, one path per series. A loading that falls towards zero is a series
leaving the index.

## What the estimator is, and is not

The model is a FAVAR whose loadings and transition coefficients both drift as
random walks:

    [ y_t ]   [ I       0     ] [ y_t ]   [  0  ]
    [ x_t ] = [ L^y_t   L^f_t ] [ f_t ] + [ e_t ],   e_t ~ N(0, V_t)

    z_t = sum_i B_it z_{t-i} + u_t,                  u_t ~ N(0, Q_t)

with `z_t = (y_t, f_t)`. The macro block is measured without error, so that part
of the state is the data and only the factors are inferred.

It is estimated by two Kalman passes — the parameters given the factors, then the
factors given the parameters — and **not** by posterior simulation. That makes it
*approximately* Bayesian, and the approximation is worth stating:

* The Kalman recursions are exact Bayesian updating, and the initial conditions
  are genuine priors.
* Nothing to do with the variances is. Both error covariances are exponentially
  weighted moving averages of past residuals, and the drift in the two
  coefficient blocks has no variance parameter at all — it is whatever the two
  forgetting factors imply. Those are plug-in point estimates substituted into
  the filter, and their uncertainty is discarded.

So the package returns one path per quantity, as `ts` objects rather than `coda`
chains, and no interval computed from them is a credible interval. Removing
exactly those priors is what makes the estimator fast — the example above runs in
about 0.01 seconds over 175 quarters and 18 series — which is what lets the
method be run over very many candidate models, and what a simulation based factor
model cannot do.

**For a fully Bayesian dynamic factor model**, with priors on the innovation
variances and posterior draws to go with them, use
[dfmtools](https://github.com/franzmohr/dfmtools) instead.

## Two deviations from the published code

* **Lags of the state.** The replication code never fills the lag blocks of the
  companion state, so it predicts from `B_1` alone. That is invisible at `p = 1`
  and wrong above it. Here the blocks are carried.
* **Missing values.** The replication code replaces them with the sample mean
  throughout. Here the filter drops a series from the periods it is missing in,
  so a series that starts halfway through the sample contributes from the date it
  starts and not before. A mean is still substituted for the principal components
  the algorithm starts from, which need a complete matrix and only set a starting
  value.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("franzmohr/fincond")
```

## References

Koop, G., & Korobilis, D. (2014). A new index of financial conditions.
*European Economic Review, 71*, 101–116.
<https://doi.org/10.1016/j.euroecorev.2014.07.002>

Raftery, A. E., Kárný, M., & Ettler, P. (2010). Online prediction under model
uncertainty via dynamic model averaging. *Technometrics, 52*(1), 52–66.

West, M., & Harrison, J. (1997). *Bayesian forecasting and dynamic models*
(2nd ed.). New York: Springer.
