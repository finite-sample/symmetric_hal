# Method and API

## Equivalent covariate codings

The primary property is independence from equivalent ways of recording a
covariate. Hold treatment and outcome definitions fixed. A recoding preserves
information when the coded value determines the original value and vice versa:
formally, a measurable bijection between supports with measurable inverse.
The conditional outcome and propensity functions compose with that inverse,
so the population ATT and its identifying assumptions are unchanged.

An estimator must state which family of such recodings it respects. Units,
zero points, and signs give coordinatewise affine maps `a + b*x`, with `b != 0`.
Arbitrary strictly monotone maps, such as logs on positive support or cubing,
preserve order up to reversal but can change all spacings. Unordered-category
permutations and invertible transformations that mix variables are separate
requirements. Binning or rounding usually loses information and is outside
this comparison.

Assigning numbers to unordered categories introduces an order that threshold
penalties can use. Even with both indicator directions and every category value
as a knot, permuting those codes can change a symmetric HAL fit. The paper's
three-category example holds the function space fixed and solves both fits
exactly, isolating the penalty's dependence on category order. The affine and
monotone guarantees do not supply arbitrary nominal-label invariance.

Full category indicators give a label-invariant construction. For one nominal
variable with declared levels, use the same schema for fitting and prediction:

```r
category_levels <- c("A", "B", "C")
category <- factor(c("B", "A", "C", "B"), levels = category_levels)
x <- model.matrix(~ 0 + category)
new_category <- factor(c("C", "A"), levels = category_levels)
newx <- model.matrix(~ 0 + category,
                     data = data.frame(category = new_category))
```

Construct the full indicator block separately for each nominal variable; a
single no-intercept formula with several factors does not generally give every
factor a full block. Relabeling categories permutes these columns. Applying the
same knot rule to every column and retaining all main effects and pairs then
permutes the functional dictionary. Both coefficient objectives and the tuned
procedure are unchanged under that permutation, with common folds and penalty
candidates. This works for both anchored and symmetric HAL. Keep declared
levels even if a training fold lacks observations in a category; new categories
outside that schema require a separate modeling decision. This construction
changes the procedure from thresholding scalar category codes. It does not
establish invariance to arbitrary contrast coding or changing an omitted
reference category. More indicator columns increase dictionary size.
See the official [R model-matrix documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/model.matrix.html)
for factor and contrast handling.

For a permitted recoding `T`, fitting the recoded training sample and predicting
at `T(newx)` should give the same result as fitting the original sample and
predicting at `newx`. Folds and other random choices must be coupled by row.
Supplied knots and variable metadata must accompany the transformation.
Agreement only at training observations is insufficient.

If held-out nuisance predictions agree, downstream clipping, calibration,
targeting, influence-function standard errors, and intervals agree whenever
they use those predictions, the unchanged outcomes and treatments, and common
settings. This is a finite-sample property. Lower prediction error is not a
prerequisite for valuing it, and invariant intervals need not have correct
coverage.

`aipw_coding_contrast(a, y, q0, g, new_q0, new_g)` in `R/att.R` compares two
effective nuisance vectors on the same observations. It returns the exact ATT
change, signed outcome and propensity components, an absolute bound, the
standard-error change and its influence-based bound, and row contributions.
The decomposition changes propensity first, then outcome predictions; it is
not a unique attribution. Supply predictions after any clipping or calibration,
with outcome predictions restored to the units of `y`. The actual calibrated
propensities determine overlap. This diagnostic applies to AIPW, not TMLE, and
equal point estimates can coexist with different standard errors.

`compare_aipw_codings(original, recoded)` in `R/causal.R` aligns saved prediction
rows and checks that folds, treatments, and outcomes match. It applies the
contrast to raw and prepared AIPW, restores outcome units, and reports both
substitution orders. `make causal-report` writes these application diagnostics
to `data/causal_coding_contrasts.csv` and generates the manuscript table.

ATT TMLE separates initial nuisance clipping from numerical bounds during
targeting. Initial outcome probabilities are clipped to `[1e-6, 1 - 1e-6]`;
targeting uses `[1e-12, 1 - 1e-12]`, with a separate propensity floor of `1e-8`.
Using the initial outcome bound during targeting can block an update and cause
the forward loss calculation to stop with a nonzero influence-function score.
The convergence flag, score residual, and step-refinement checks must all pass.
`q_targeting_floor` in `prepare_att()` and `compare_att()` controls the numerical
outcome bound; it must be positive, numerically distinguishable from one after
subtraction, and smaller than the initial `q_bound`.

Requiring pointwise invariance under every bijection would be too broad for
using numerical order or distance: swapping two unseen values while fixing the
training values would force their predictions to agree. The choice of admissible
recodings therefore needs a substantive interpretation.

## Symmetric zero-order HAL

Symmetric HAL fits an intercept plus zero-order indicator products. Each
indicator can point up (`x >= knot`) or down (`x <= knot`); products involve at
most two distinct covariates. The intercept is unpenalized and features are
not standardized.

For Gaussian outcomes the objective is mean squared residual divided by two,
plus `lambda * sum(abs(beta))`. For binomial outcomes, replace squared error
with mean negative log likelihood. The penalty is strictly positive.

## A uniquely defined fitted function

A lasso can have several coefficient vectors with the same training predictions
and objective. Those vectors need not agree at new covariates. The estimator
selects the vector of minimum Euclidean norm among all primary minimizers.
The squared norm is strictly convex, so this secondary vector is unique.

Reflecting coordinates permutes the full symmetric functional dictionary.
Both objectives are unchanged by that permutation, which carries the unique
selected solution to its reflected version. Predictions at correspondingly
reflected inputs agree. The same argument covers cross-validation when folds,
penalty candidates, and the tuning rule are shared.

This argument requires retaining functions that agree on the observed sample
but differ elsewhere. Sample-column deduplication is therefore inappropriate.
The default knot rule uses type-7 quantiles at probabilities `k/(q + 1)` for
`k = 1, ..., q`; the probability grid is symmetric. Duplicate knot *values*
within a coordinate are removed because they define identical functions.
Equality is included in both indicator orientations. An explicit `knots` list
can instead supply fixed numeric thresholds for each coordinate. Transform
those thresholds along with the covariates when comparing codings.

The same proof covers all coordinatewise affine recodings, not just signs:
type-7 quantiles commute with positive affine maps and reverse quantile indices
under negative scales. This guarantee assumes exact arithmetic and optimization;
finite-precision behavior is checked separately.

It does not cover arbitrary nonlinear monotone maps. With training values
`c(0, 2)`, three interpolated knots are `c(.5, 1, 1.5)`. After cubing the
values, the knots become `c(2, 4, 6)`, rather than the cubes of the original
knots. The training design stays identical, but prediction boundaries move.
The paper gives an exact lasso example, without simulation or solver error.

For the stronger monotone guarantee, set `knot_rule = "observed"` with
`symmetric = TRUE`. For each training size `n` and requested count `q`, select
ranks `r_k = 1 + floor((n - 1)*k/(q + 1))` for `k = 1, ..., q`, together with
every reversed rank `n + 1 - r_k`. Use the corresponding observed values,
removing repeated values. This rule is independent of spacings and selects
at most `2*q` distinct knots per coordinate. For an even sample and `q = 1`,
it retains both middle observations rather than their average. Ties and
constant coordinates are allowed. The rule is rebuilt in every validation
training fold and in the final fit.

The paper proves functional equivariance of the ideal symmetric estimator
under arbitrary coordinatewise strictly monotone recodings with this rule.
Such a predictor is constant within gaps between observed values. Gaussian
and binomial regression tests check nonlinear increasing and decreasing
recodings, validation predictions, tuning, and new prediction points.
Numerical optimization retains the same certificate checks as the default.
The empirical studies use the quantile default; their accuracy comparisons
do not evaluate the observed-rank option. Supplied `knots` bypass the default
rule and cannot be combined with `knot_rule = "observed"`.

Full pointwise monotone invariance also restricts continuity. On `R^d`, or a
product of continuous intervals, a fixed-sample predictor continuous in every
coordinate and invariant to all increasing recodings must be constant in the
query. Continuity connects the constant predictions on either side of every
observed value. Thus a nonconstant continuous spline fit cannot meet this full
requirement, although it can respect affine recodings. The observed-rank
indicator construction permits jumps. This restriction does not apply in the
same way to discrete or disconnected support, and says nothing against smooth
estimation or population identification generally.

## Fitting

```r
fit_hal(x, y, lambda, num_knots = 3, symmetric = TRUE,
        family = "gaussian", solver_tol = 1e-20,
        certificate_tol = 1e-6, face_tol = 1e-8, knots = NULL,
        knot_rule = "quantile")
```

1. Build the full dictionary and solve the primary lasso using `glmnet`.
   Binomial fits follow a geometric path from the intercept-only penalty to
   the requested penalty. Only the requested endpoint enters the estimator;
   the intermediate penalties are computational warm starts.
2. Identify its equicorrelation/sign face, allowing for measured primary
   optimization error.
3. Use `quadprog` to minimize coefficient norm while retaining the primary
   fitted linear predictors and the face's sign restrictions.
4. Check primary and secondary stationarity, feasibility, and preservation of
   the primary predictions, L1 norm, and objective. The secondary check includes
   every dictionary coordinate, including those excluded from the proposed face.

The returned list contains `beta`, `intercept`, `lambda`, `features`, `fitted`,
`objective`, `primary_objective`, `diagnostics`, `family`, and `predict`.
Use `fit$predict(newx)` for response predictions or
`fit$predict(newx, type = "link")` for linear predictors.

A convergence failure or residual above `certificate_tol` stops the fit.
The tolerance is absolute. It is a numerical optimality criterion, not a
universal bound on out-of-sample prediction error. Tests therefore compare
held-out predictions separately, including data with ties, duplicate
covariates, and different outcome scales.

## Cross-validation

```r
cv_hal(x, y, lambda, foldid, num_knots = 3, symmetric = TRUE,
       family = "gaussian", tie_tol = 1e-8,
       solver_tol = 1e-20, certificate_tol = 1e-6, knots = NULL,
       knot_rule = "quantile")
```

`lambda` is an explicit vector of positive candidates. `foldid` assigns every
row to one of at least two integer-labelled folds. Each training fold rebuilds
its own dictionary; validation outcomes do not determine knots or fitting.
Binomial training sets must contain both 0 and 1.

Candidate loss is averaged over validation observations, so uneven fold sizes
are weighted by their sizes. Among penalties whose loss is within `tie_tol`
of the minimum, select the largest, then refit on all rows. The returned fit
also contains `cv` (candidate penalties and losses), `cv_diagnostics` (every
fold/candidate's residuals), `foldid`, out-of-fold link predictions in
`validation_predictions`, and validation/refit times in `seconds`.

Finite-precision losses near the tuning threshold can still select different
candidates under perturbation. Reusing folds and checking the resulting
predictions are necessary numerical checks alongside the exact theorem.

## Canonical orientation

Source `R/canonical.R` after `R/basis.R` and `R/fit.R` to use
`cv_canonical_hal(x, y, lambda, foldid)`. It uses the anchored dictionary after
orienting each covariate by the sign of its normalized centered third moment.
Coordinates with absolute score at most `1e-12` retain both directions; their
fits are averaged. Binary predictions are averaged as probabilities.

The rule is recomputed inside each validation training fold. It uses one
orientation on generic samples, but up to `2^d` on fully ambiguous samples.
The returned fit includes `signs`, constituent `fits`, per-fold `cv_signs`,
and the usual validation losses, predictions, numerical diagnostics, and times.
This gives a simpler invariant comparator; it does not combine atoms from
different orientations in one penalized fit. Exact-arithmetic equivariance
does not remove the need to check numerical behavior near orientation or
validation-loss thresholds.

The third-moment rule uses spacings and is not invariant under arbitrary
increasing recodings. More generally, an unsupervised column-only orientation
rule that is row-permutation invariant cannot select one direction on distinct
observations while respecting all increasing recodings and sign reversals.
An increasing map can send the sorted values to the reversed negatives; the
result is the reflected sample up to row permutation. Such a rule must retain
both signs. Applying it separately in each coordinate would require all
`2^d` orientations. This restriction does not cover joint or outcome-informed
orientation rules.

For the first-order package comparison, `fit_nuisance("hal9001_canonical", ...)`
standardizes and orients each outer-training nuisance sample before invoking
`hal9001` with its default basis and shared inner folds. Ambiguous orientations
are tuned separately by the package, then their response predictions are
averaged. Thus this wrapper preserves the package's fitting procedure after
preprocessing; it differs from the fold-local, common-penalty zero-order
comparator above. Both procedures keep outer held-out observations out of
nuisance training. Neither supplies a guarantee of causal accuracy.

## Computational limits

The dense design and quadratic program make this implementation suitable for
small and moderate dictionaries. No first-order basis, empirical screening,
or relaxed refit is implemented. Absolute tolerances and conditioning can
cause valid mathematical problems to fail numerically. Do not suppress such
failures in comparative studies or substitute an arbitrary lasso solution.

The anchored option uses the same optimization procedure with positive
orientations only. This provides a controlled comparison of dictionaries;
`hal9001` defaults are a separate software baseline with different modeling
and tuning choices.

## Known-effect study reports

Collect complete run directories, then generate reports without refitting
nuisance models or drawing new samples:

```sh
Rscript R/semisynth_collect.R results_completed_main results_completed_canonical
Rscript R/semisynth_report.R
```

Replace the directory names with the completed runs. Collection recomputes
targeting from saved nuisance predictions using the current procedure. The run
stage saves nuisance predictions; collection estimates the treatment effect.
Set `SEMISYNTH_CORES` to parallelize that step.
Collection checks the
frozen population model, sample size, recodings, paired observations and folds,
and duplicate or failed jobs. Reporting requires every method and estimator
on the same replicate grid under both codings. It produces:

- `semisynth_summary.csv`: bias, RMSE, Monte Carlo standard errors, and exact
  binomial intervals for empirical coverage.
- `semisynth_comparisons.csv`: paired MSE differences and RMSE ratios for the
  planned method comparisons, with pointwise approximate Monte Carlo intervals.
- `semisynth_movement.csv`: coding movement in ATT estimates and standard errors.
- `semisynth_nuisance_summary.csv`: nuisance prediction losses and drift summaries.
- `semisynth_diagnostics.csv`: preparation clipping, nuisance warnings, and times,
  counted once per fit rather than once per reported estimator. Times distinguish
  nuisance fitting from the newly computed targeting step.

All files are written to `data/`. The `prepared_clipped_*` columns refer to the
preparation shared by calibrated AIPW and TMLE; raw AIPW uses different inputs.
The ratio intervals use a paired delta calculation on the log scale and a
`t` critical value. Log-ratio intervals are undefined when either empirical MSE
is zero; the point ratio is zero for a zero numerator and positive denominator.
Twenty exploratory draws cannot establish interval calibration or a general
accuracy ranking. The method-comparison intervals are not intervals for the
causal effect, and recorded wall times reflect concurrent machine load.

The nuisance losses use the known functions at the sampled covariate rows.
`q0_risk` and `g_risk` are mean squared errors of the unprepared predictions.
`q0_treated_risk` and `q1_treated_risk` weight squared errors by the true
treatment probabilities, normalized to sum to one. They assess the treated
covariate distribution rather than restricting evaluation to the observations
whose realized treatment is one.

For effective nuisance predictions `q` and `g`, the drift diagnostic is

```r
mean((g - true_g) * (q - true_q0) / (1 - g)) / mean(true_g)
```

`raw_drift` uses the original outcome predictions and propensities capped at
0.99. `prepared_drift` uses the clipped outcome predictions in outcome units
and the calibrated propensities. For fixed nuisance functions on a fixed
covariate distribution, this expression equals the mean ATT estimating score
at the true effect. The reported quantity substitutes the sampled covariate
distribution and the fitted predictions. It is a diagnostic of the joint
nuisance error and overlap, not the finite-sample bias of the fitted estimator:
the fits and preparation depend on the data, and the estimator uses the random
observed treated count. It does not measure TMLE bias after targeting.
