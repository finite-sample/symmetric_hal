# The Highly Adaptive Lasso depends on how you code your covariates

The zero-order Highly Adaptive Lasso (HAL) regresses the outcome on indicators
`1{x_j >= u}` and their products, with an L1 penalty. Every indicator switches
on as a covariate *increases*, so a function that switches on as a covariate
decreases has to be built from differences:

```
1{x1 <= a, x2 <= b} = 1 - 1{x1 > a} - 1{x2 > b} + 1{x1 > a, x2 > b}
```

That costs three units of penalty; `1{x1 >= a, x2 >= b}` costs one. So HAL's fit
changes when you code a covariate as age instead of years to retirement, or
income instead of a poverty indicator. Where HAL is the nuisance estimator in
TMLE or double machine learning, the treatment-effect estimate moves too.

**Symmetric HAL** adds the reflected indicators `1{x_j <= u}`. Its fit is
exactly the same under any sign recoding of any covariate, and it keeps HAL's
`n^(-2/3)` rate.

## Use it with hal9001

No new package is needed. Give `hal9001` the covariates twice, as `x` and `-x`;
its indicator `1{-x_j >= -u}` is `1{x_j <= u}`. `R/basis.R` builds the
symmetric basis and translates it into a `hal9001` basis list:

```r
source("R/basis.R")

fit_shal <- function(x, y) {
  f <- make_dictionary(x, num_knots = 9, symmetric = TRUE)
  hal9001::fit_hal(cbind(x, -x), y, basis_list = to_hal9001(f, ncol(x)),
                   smoothness_orders = 0, max_degree = 2)
}
fit <- fit_shal(x, y)
predict(fit, new_data = cbind(newx, -newx))
```

All of `fit_hal`'s families and cross-validation options work unchanged. The
test suite checks that this route gives the same fitted values as a direct
lasso on the symmetric dictionary.

## What we find

On real data, negating covariates at random (19 distinct codings, same folds
and test split):

| | HAL | symmetric HAL |
|---|---|---|
| LaLonde (185 NSW trainees vs PSID-1): effect of training on the trainees' 1978 earnings, `hal9001` default first-order fits | -$16,526 to $1,849 (17.8 standard errors) | $1,852 to $1,903 |
| Same, zero-order fits | $451 to $3,551 (2.6 standard errors) | $2,800 |
| LaLonde, default first-order fits, earnings in dollars vs thousands and age in years vs months | $1,181 vs $47,975 | $1,639 in both (scaled) |
| Boston housing: typical change in a test prediction | 0.11 SD of the outcome | under 0.001 SD |
| Boston housing: test RMSE | 2.59 to 3.26 | 2.97 |

In a semi-synthetic LaLonde with a known effect of $1,800 (66 replicates),
the default `hal9001` estimator under the original coding is biased by -$2,058
and its 95% intervals cover 36% of the time; negating the continuous covariates
removes the bias. Symmetric first-order HAL is unbiased in both codings (-$17,
-$19), at 0.95 times the efficiency bound, with 97% coverage.

The randomized experiment puts the LaLonde effect at $1,794 (SE $671). Agreement
with one experiment says little about accuracy, and symmetric HAL is not more
accurate than a typical coding of HAL on Boston. What it removes is the
dependence on the coding.

In simulations (n = 500, 10 covariates, 50 replicates per cell), symmetric HAL
beats HAL by up to 2.3x in risk when a box points away from HAL's anchor. When
the box points toward it, HAL is 3 to 18% better, the price of a dictionary
four times larger for interactions.

HAL is not alone. Refitting 18 common regression estimators on Boston housing
under recodings that carry the same information (negation, shifts, units, the
reference level of a categorical covariate), least squares, the additive model
and CART never move. MARS (`earth`), histogram-based boosting (`xgboost`), the
lasso with interactions (under shifts), first-order HAL (under units), and every
penalized or distance-based estimator given dummy variables (under the reference
level) do. Symmetric HAL removes the one-sided-basis problem; with first-order
hinges, also divide each covariate by its standard deviation, or the fit depends
on units.

## Reproduce

```
make check    # lint and tests
make report   # figures and tables from the collected results in data/
make paper    # report, then compile paper/ms.pdf
make sims     # rerun the simulations (hours; CORES=8 REPS=50)
make recode   # rerun the prediction recodings
make lalonde  # rerun the LaLonde recodings (hours)
make semisynth  # rerun the semi-synthetic LaLonde simulation (hours)
make audit    # rerun the 18-estimator coding audit
make data     # re-download the Dehejia-Wahba and PSID-1 files
```

R dependencies are listed in `DESCRIPTION`
(`remotes::install_deps(dependencies = TRUE)`).

| Path | Contents |
|---|---|
| `R/basis.R` | anchored and symmetric dictionaries; `to_hal9001()` |
| `R/fit.R` | the fitters used in the simulations |
| `R/sim_*.R` | simulation design, run, collect, report |
| `R/nuisance.R`, `R/att.R` | HAL nuisance fits (zero order and `hal9001`, anchored and symmetric); cross-fitted effect on the treated |
| `R/recode_*.R` | prediction recodings: run, collect, report |
| `R/lalonde_*.R` | LaLonde data, recodings: fetch, run, collect, report |
| `R/semisynth_*.R` | semi-synthetic LaLonde: run, collect, report |
| `R/audit_*.R` | coding audit of 18 estimators: run, collect, report |
| `data/` | the LaLonde data and collected results, one row per fit |
| `paper/` | manuscript and generated tables |
| `figures/` | generated figures |
| `tests/` | checks of every proposition and of the hal9001 route |
