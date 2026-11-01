# Covariate coding and causal estimates with HAL

Changing covariate signs can change a treatment-effect estimate even though
the data contain the same information. This project reproduces that sensitivity
for a specified HAL nuisance-estimation pipeline and compares practical remedies.
The [causal estimates](data/causal_summary.csv) document the observational
finding. The [AIPW coding contrasts](data/causal_coding_contrasts.csv) connect
changes in saved nuisance predictions to changes in estimates and standard
errors. A completed [known-effect causal pilot](data/semisynth_movement.csv)
checks coding movement on repeated samples. Prediction and causal accuracy
comparisons assess separate tradeoffs.

The code includes symmetric HAL and a simpler canonical-orientation comparator.
Symmetric HAL includes both directions of each indicator and selects the
coefficient vector of smallest Euclidean norm among lasso solutions. Canonical
orientation chooses directions from training covariates and averages ambiguous
directions. Both support Gaussian and binomial outcomes, main effects, and
pairwise interactions. Invariance alone establishes no accuracy advantage.
The primary benefit is a result independent of the permitted coding choices;
it does not require lower prediction risk. The ideal symmetric estimator
is equivariant to coordinatewise changes of units, zero points, and signs with
exact arithmetic and exact optimization. Numerical behavior is checked
separately. Arbitrary monotone recodings require a stronger condition,
which interpolated quantile knots do not meet. The explicit
`knot_rule = "observed"` option uses observed ranks to meet that stronger
requirement for the ideal symmetric estimator. The existing studies use the
quantile default. The
[method guide](docs/method.md) states this distinction explicitly.

The [paper](paper/ms.tex) defines equivalent codings and the required agreement
of predictions and inference, then develops the causal example, remedies,
and penalty geometry. Risk bounds are collected in an appendix.
The [method guide](docs/method.md) explains
the fitting procedure. The [research plan](docs/research_plan.md) describes the
mechanism-based checks and study designs.

## Get started

Use R 4.4 or later with glmnet 5.0 or later. From this directory:

```r
install.packages(c("Matrix", "glmnet", "quadprog"))
source("R/basis.R")
source("R/fit.R")

set.seed(7)
x <- matrix(runif(160, -1, 1), ncol = 2)
y <- 2 * (x[, 1] <= -0.2 & x[, 2] <= -0.2) + rnorm(80, sd = 0.6)
foldid <- sample(rep_len(1:5, nrow(x)))

fit <- cv_hal(x, y, lambda = c(0.01, 0.03, 0.1), foldid = foldid)
newx <- rbind(c(-0.8, -0.8), c(0.8, 0.8))
fit$predict(newx)
fit$lambda
fit$diagnostics
```

For a fixed penalty, use `fit_hal(x, y, lambda = 0.03)`. Set
`symmetric = FALSE` for the anchored dictionary with the same solution-selection
rule. For binary outcomes, set `family = "binomial"`; predictions are
probabilities, or use `fit$predict(newx, type = "link")` for log odds.
Every binomial training fold must contain both outcomes.

Both fitting functions default to three interior quantile knots per coordinate.
`x` and prediction inputs must be finite numeric matrices with columns in the
same order. `y` must be a finite numeric vector. Missing data and categorical
encoding must be handled before fitting.
The invariance guarantee does not cover arbitrary permutations of numeric codes
for unordered categories; the paper gives an exact three-category example.
For nominal variables, use one indicator per category and keep the full category
schema fixed across training folds and prediction inputs. Relabeling categories
then only permutes columns, which preserves the ideal fitted procedure when
the knot rules, tuning candidates, and folds agree. This applies to both
anchored and symmetric dictionaries. Omitting a reference category or changing
contrast bases does not have this guarantee. The
[method guide](docs/method.md) gives an encoding example.

For invariance under arbitrary coordinatewise strictly monotone recodings,
use `knot_rule = "observed"` in `fit_hal()` or `cv_hal()` with the symmetric
dictionary. Knots are observed order statistics, selected symmetrically in
rank and rebuilt inside each training fold. This can use up to twice
`num_knots` distinct thresholds per coordinate. The
[method guide](docs/method.md) gives the rank rule and its exact-arithmetic
guarantee; regression tests check held-out prediction agreement.

## Compare sign codings

Reuse the penalty candidates and fold assignments, and transform prediction
inputs in the same way as training inputs:

```r
signs <- c(-1, 1)
reflected <- cv_hal(sweep(x, 2, signs, `*`), y,
                    lambda = c(0.01, 0.03, 0.1), foldid = foldid)
max(abs(fit$predict(newx) -
        reflected$predict(sweep(newx, 2, signs, `*`))))
```

In exact arithmetic, the symmetric estimator gives the same predictions.
Numerical fits must pass checks of both optimization stages, and the test suite
checks held-out predictions under reflection, including tied and duplicate
covariates. A failed convergence or optimality check raises an error.

For the canonical-orientation comparator, using the same data and folds:

```r
source("R/canonical.R")
canonical <- cv_canonical_hal(x, y, c(0.01, 0.03, 0.1), foldid)
canonical$predict(newx)
canonical$signs
```

The orientation rule is recomputed within each training fold. Generic samples
use one orientation; ambiguous coordinates require averaging additional fits.

## Reproduce the examples and paper

Install the project dependencies listed in `DESCRIPTION`. Install pak with
`install.packages("pak")`, then run
`pak::local_install_deps(dependencies = TRUE)`. The paper also needs a LaTeX
installation with `latexmk`.

To refit the reported studies, use R 4.6.0 and the core analysis versions below
after installing the project dependencies. In a fresh R session:

```r
dir.create(".Rlib", showWarnings = FALSE)
.libPaths(c(normalizePath(".Rlib"), .libPaths()))
install.packages("remotes", lib = ".Rlib", repos = "https://cloud.r-project.org")
reference <- c(Matrix = "1.7-5", glmnet = "5.0", quadprog = "1.5-8",
               hal9001 = "0.4.6", tmle = "2.1.1")
for (package in names(reference)) {
  remotes::install_version(package, version = reference[[package]],
                           lib = ".Rlib", upgrade = "never",
                           repos = "https://cloud.r-project.org")
}
stopifnot(all(vapply(names(reference), function(package) {
  packageVersion(package) == package_version(reference[[package]])
}, logical(1))))
```

The standard Make targets use `.Rlib` first. For direct `Rscript` commands,
set `R_LIBS=.Rlib`. The
[version installer](https://remotes.r-lib.org/reference/install_version.html)
can retrieve archived source packages; compiling them requires R's development
tools. These pins restore the core fitting and targeting packages, not the
entire operating environment. The committed
[nuisance-run](data/causal_session.txt) and
[targeting-session](data/causal_targeting_session.txt) records give the reference
R, platform, and dependency versions. Numerical comparisons remain subject to
the stated tolerances across platforms.

```sh
make check
make example
make paper
```

These commands run lint and scientific regression tests, regenerate the small
examples, and build the paper, respectively. The paper build regenerates its
reports from committed data without refitting the long simulation studies.
The examples use fixed seeds and explicit tuning grids. The causal section uses
committed out-of-fold predictions to reproduce ATT estimates under three covariate
codings. Neither example establishes superior prediction risk or causal accuracy.

## Scope

This is a research implementation for small and moderate dictionaries. The
secondary optimization uses dense matrices. With `d` covariates and `q` distinct
knots per covariate, the symmetric dictionary has
`2*d*q + 4*choose(d, 2)*q^2` columns. Start with a small knot count.

The code retains distinct basis functions even if their observed columns agree.
It does not provide first-order splines or screening. Minimum-norm selection
resolves lasso nonuniqueness; adding reflected columns alone is insufficient for
invariance away from the training sample. Numerical tolerances are absolute and
can require attention when the outcome has a large scale.

The theoretical risk bound concerns a constrained estimator on a full lattice
with fixed dimension and Gaussian errors. It is not a risk guarantee for this
cross-validated finite-knot implementation. The prediction pilot below compares
mixed orientations with selecting or averaging anchors and finds no general
risk advantage.

## Treatment-effect coding sensitivity

In the LaLonde/PSID example, reversing continuous covariate signs changes the sign
of an ATT TMLE using first-order `hal9001` nuisance fits. The sample, folds, and
estimation settings stay fixed. The [generated estimates](data/causal_summary.csv)
include AIPW and TMLE; [targeting diagnostics](data/causal_estimates.csv) record
convergence, score residuals, and step-size sensitivity. The different-units case
produces nuisance convergence warnings, recorded in
[fit diagnostics](data/causal_fits.csv).

This is a finding about coding sensitivity of a specified pipeline. It does not
identify which estimate is correct, establish interval coverage, or demonstrate
that symmetric HAL improves causal accuracy. The causal application uses the
package's first-order defaults; the symmetric implementation here is zero-order.

```sh
make causal-report
Rscript R/causal_run.R results_new_causal 3
Rscript R/causal_collect.R results_new_causal
make causal-report
```

The first command recomputes estimates from committed nuisance predictions. The
remaining commands refit nuisances and replace those predictions. Use a fresh run
directory when source code or dependencies change. The
[nuisance-fitting session](data/causal_session.txt) records the environment used
to fit the saved predictions. The
[targeting session](data/causal_targeting_session.txt) records the environment
and source hashes for recomputing estimates from those predictions.

## Research experiments

`R/pilot_run.R` compares anchored, selected-anchor, averaged-anchor, and symmetric
HAL, with CART and random forests as secondary benchmarks. It records all planned
cells, fitting failures, numerical certificates, and paired prediction risks. The
[research plan](docs/research_plan.md) specifies the grid and decision criteria.
The completed prediction pilot finds mixed risk differences, with no general
advantage for symmetric HAL. In the fixed-knot comparison, its risk relative to
selected anchors changes direction with sample size. Large gains in some discrete
cases can also reflect a larger dictionary span, not just a different penalty.
The paper reports the paired Monte Carlo uncertainty and all planned scenarios.

To refit the complete prediction pilot, including canonical orientation, and
regenerate its reports:

```sh
make pilot CORES=4
```

This is a long computation. `PILOT_DIR` and `CANONICAL_DIR` can select fresh
checkpoint directories. `R/canonical_run.R` can also add the comparator to an
existing primary run; its optional final argument selects job IDs for feasibility.
The collector requires both runs to be complete and checks matching folds and
outcome scales before producing paired comparisons.
After collecting the complete 20-replicate studies, `Rscript R/pilot_report.R`
generates one risk-comparison figure per signal and the fixed-knot results table.
Its intervals describe Monte Carlo uncertainty in paired risk ratios.
Fitting times are elapsed times under concurrent machine load, not isolated
serial benchmarks. `data/pilot_diagnostics.csv` records the observed prediction
movement under reflection and the canonical orientation counts.

`Rscript R/causal_path_check.R results_causal_controls` reproduces a numerical
diagnostic in the tighter first-order control. One inner validation path ends
before the selected master penalty and reuses its endpoint prediction. Extending
that path changes validation predictions even though the selected master fit
is reproducible. This diagnoses that control; it does not establish the cause
of the default-coding reversal.


`R/semisynth.R` defines a known-effect population using the empirical covariate
distribution and Gaussian outcomes calibrated to the earnings data. It fixes
the ATT at 1,800 and saves the true nuisance functions. `make semisynth` compares
AIPW and ATT TMLE under original and reflected codings, including oracle-initialized
estimates. The completed 20-draw exploratory pilot finds exact coding agreement
for both canonical procedures and agreement to numerical precision for symmetric
HAL. Its [paired accuracy comparisons](data/semisynth_comparisons.csv) vary by
causal estimator and do not establish a general accuracy advantage. Gaussian
outcomes can be negative;
the simulation does not claim to reproduce the full earnings distribution.
The `hal9001_canonical` nuisance method also applies standardized canonical
orientation directly to the first-order package. It uses outer-training
preprocessing followed by package cross-validation, with separately tuned fits
averaged when directions are ambiguous. The method guide gives the distinction
from the zero-order comparator. The semi-synthetic runner includes it by default;
pass it explicitly to the causal runner to reproduce the application comparison.
In the observational application, it gives the same ATT TMLE across the three
codings to numerical precision, while retaining nuisance convergence warnings.
The completed solver and scaling controls preserve the point-estimate sign
reversal. Standardization removes the observed unit dependence to numerical
precision at both solver settings; tightening the solver alone does not.
Nuisance convergence warnings remain, so these comparisons do not rule out
optimization effects. `make causal-control-report` regenerates both comparison
tables from the committed `data/causal_control_estimates.csv`.

To refit the solver, scaling, and canonical-orientation controls:

```sh
Rscript R/causal_run.R results_new_controls 3 hal9001_tight,hal9001_scaled,hal9001_scaled_tight
Rscript R/causal_run.R results_new_canonical 3 hal9001_canonical
Rscript R/causal_control_report.R results_new_controls results_new_canonical
```

The report requires both completed run directories and recomputes targeting
from their saved nuisance predictions. Both fresh collection and reporting
from cached estimates require all three numerical controls and the canonical
comparison under all three codings. An incomplete comparison stops the build
before it writes tables. Nuisance warnings are reported separately from
targeting checks.

`make semisynth-paper` generates the completed causal pilot's accuracy, coding-movement,
and paired-comparison tables from the collected estimates. It recomputes
summaries and rejects an incomplete method or replicate set. `make paper` includes
this target. The tables keep
Monte Carlo uncertainty separate from uncertainty about the treatment effect.
Collection also writes `data/semisynth_session.txt`, recording the nuisance-run
source fingerprints and environments alongside the current targeting sources
and environment.

`Rscript R/semisynth_oracle.R` reproduces the
[oracle benchmark](data/semisynth_oracle_summary.csv). It compares empirical
influence-function intervals with intervals using the known conditional variance
over independent simulation draws. That conditional variance is available only
for this raw oracle AIPW calculation, not for fitted nuisance or TMLE estimates.

`R/sim_setup.R` supplies simulation truths and paired fitting functions.
`R/recode_run.R` compares codings on real regression data. The LaLonde and
semi-synthetic runners study nuisance estimation; `data/lalonde_dw.csv` contains
the empirical input. These are research drivers, not evidence of efficacy.

The long-run targets `make sims`, `make recode`, `make lalonde`,
`make semisynth`, and `make audit` write private checkpoints under
`results_current_*/`. Run a small feasibility study before the full targets;
the dense solver can become costly. Each run records its complete job grid, source-input fingerprints, package versions, and R
session. Collectors reject missing or failed jobs. Use a fresh checkpoint
directory after changing code, configuration, or the R environment.
`make benchmark-reports` requires completed outputs from all five studies;
it is separate from the paper build. Follow the research plan before launching
comparative studies or making performance claims.

## Files

- `R/basis.R`, `R/fit.R`: dictionary, fitting, and cross-validation.
- `R/canonical.R`, `R/nuisance.R`: canonical orientation and nuisance learners.
- `tests/testthat/`: optimization, equivariance, and penalty-geometry checks.
- `R/example.R`, `R/fig_coding.R`: reproducible examples.
- `paper/`: manuscript and generated table fragments.
- `docs/`: method details and research design.

MIT license. See [LICENSE](LICENSE).
