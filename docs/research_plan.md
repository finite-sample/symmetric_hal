# Research plan: coding sensitivity and practical remedies

The primary requirement is that equivalent expressions of the covariates give
one analysis: predictions, treatment-effect estimates, standard errors, and
intervals should agree under a declared family of information-preserving
recodings. Lower prediction risk is a separate property, not an admission
criterion for an invariant method.

## Formalize before choosing experiments

1. Hold outcome and treatment definitions fixed. Define a recoding as a
   measurable bijection between covariate supports with measurable inverse.
   Show that the conditional nuisance functions compose with that inverse and
   the ATT functional is unchanged. Exclude information-losing coarsening.
2. Declare the admissible transformations. Distinguish coordinatewise affine
   unit/origin/sign changes from arbitrary strictly monotone recodings,
   nominal-category permutations, and joint transformations such as rotations.
   Information preservation alone cannot select the family: unrestricted
   pointwise bijection invariance forces equal predictions at unseen values.
3. Require equality at correspondingly transformed prediction points for every
   training sample, with row-attached folds and auxiliary randomness coupled.
   Prove the transfer from invariant held-out nuisances to the downstream
   estimate and interval. Equality of a scalar ATT alone is a weaker property.
4. Establish the exact mathematical guarantee before numerical assessment.
   For symmetric HAL, verify functional dictionary closure, penalty invariance,
   solution selection, fold-local knots, and tuning ties. Separate that ideal
   estimator from finite-precision implementation.
5. Derive targeted checks from each possible failure: changed function space,
   changed penalty on a common function space, nonunique off-sample solutions,
   tuning differences, or incomplete numerical optimization. Use deterministic
   counterexamples and controlled algebra before population simulations.

The ideal symmetric estimator defined by the current dictionary and exact
two-stage optimization meets the affine requirement in exact arithmetic:
its symmetric type-7 probability grid transforms correctly under
positive and negative affine maps. It does not meet the full monotone
requirement. The paper's cubic-recoding example has identical training design
matrices but different held-out predictions because interpolated knots move
within gaps. Fixing rank indices by training size, independently of observed
spacings, and closing them under rank reversal would remove that source of
dependence when the corresponding order statistics are used as knots.
The manuscript states this construction as a corollary, including ties,
constant coordinates, and fold-specific training sizes. Its necessary-condition
proposition also shows that full increasing-recoding invariance forces
predictions to be constant between observed coordinate values, including the
exterior intervals. That restriction alone is insufficient: the heights of
the pieces must not depend on spacings either. Deterministic checks use supplied
rank knots to compare complete dictionaries and held-out Gaussian and binomial
predictions under increasing and decreasing nonlinear recodings. The explicit
`knot_rule = "observed"` option implements the construction, including rank
selection inside each validation fold. End-to-end Gaussian and binomial checks
compare tuning losses, selected penalties, validation predictions, and new
predictions under nonlinear recodings. The public default and all existing
study specifications remain quantile-based; their risk results do not evaluate
the observed-rank option.

On `R^d`, or a product of continuous intervals, the gap restriction plus
continuity forces a fixed-sample predictor to be constant in the query. Thus
full monotone invariance and nonconstant continuous spline prediction cannot
both be required. This is a scope condition for the stronger guarantee, not
a criticism of smooth estimation or a consequence of information loss.
Distinguish the affine requirement from the stronger order-based requirement
when interpreting comparisons between spline and indicator learners.

The stronger requirement also changes the canonical comparator. A column-only,
row-permutation-invariant orientation rule that ignores increasing recodings
must retain both directions on distinct observations if it is to respect
reversal. Its Cartesian product then has `2^d` orientations. The third-moment
rule avoids this restriction by using spacings. Compare this marginal-frame
cardinality with the symmetric dictionary's polynomial number of columns at
fixed interaction degree, without claiming a universal computational lower
bound: joint or outcome-informed rules and shared fitting can behave differently.

## Experiments implied by the formalization

Nominal categories require a separate construction. The exact three-level
example shows that permuting scalar category codes changes the threshold
penalty even with a saturated function space. Full category indicators turn
relabeling into a column permutation and preserve the ideal anchored or
symmetric procedure with identical column rules, folds, and tuning candidates.
Keep the declared schema fixed, including levels absent from a training fold.
Deterministic Gaussian and binomial checks verify the tuned construction;
arbitrary contrasts and omitted-reference changes are outside this guarantee.

For procedures satisfying the theorem, coding movement should be zero for any
sample size, truth, or signal-to-noise ratio in exact arithmetic. Numerical
checks should therefore stress ties, duplicate columns, near-ambiguous canonical
scores, large or small units, and nearly tied tuning losses. These checks assess
implementation fidelity, not whether a mathematical invariance holds on average.

For noninvariant procedures, compare the same sample under permitted recodings.
First hold the dictionary span, candidate penalties, folds, and solution rule
fixed wherever the mechanism allows. Penalty examples identify functions whose
regularization changes with anchor; knot examples identify changes in the
available functions; solver controls isolate incomplete paths. Do not interpret
one mechanism's experiment as identifying another.

Population simulations can then quantify the frequency and magnitude of coding
movement, and assess the accuracy and computation of invariant alternatives.
Vary sample size, regularization, signal geometry, overlap, or knot density only
with a stated mechanism and a paired design. The existing broad prediction
pilot and completed causal accuracy pilot remain exploratory evidence. Do not
expand them before settling the transformation family and the question each
additional axis answers.

For AIPW, the paper derives an exact coding contrast from the effective nuisance
vectors. Outcome-prediction changes are weighted by treatment and control odds;
propensity changes are weighted by control residuals and the inverse product
of the two control probabilities. This identifies overlap, residual size, and
the treated fraction as possible amplifiers. Check signed contributions and
their cancellation before treating a worst-case bound as a typical effect.
Compare influence values and standard errors as well as point estimates:
the paper gives an exact example with identical estimates but different standard
errors. Apply the identity to prepared predictions in outcome units, using the
actual calibrated propensities. Initial clipping need not bound them after
calibration. These AIPW identities do not explain TMLE movement without further
analysis of its targeting map.

The first-order spline calculation identifies another controlled mechanism.
With a free intercept and penalized linear term, the univariate anchored
penalty is the absolute slope at the anchor plus the total variation of the
derivative. Reflection changes the endpoint slope term while preserving the
function space. The paper gives exact, distinct lasso solutions for a
three-observation sample and verifies the package's corresponding basis atoms.
This separates penalty dependence from numerical optimization and does not
attribute the observational TMLE reversal to that mechanism. Leaving affine
functions unpenalized removes this particular univariate penalty difference;
it is not a general multivariate remedy. Any subsequent controlled comparison
should specify whether boundary terms, knots, tuning, or optimization change.

## Secondary prediction comparison

Add a canonical-orientation comparator on the identical prediction draws and
folds. Within each validation training fold, center each covariate, divide by
its maximum absolute centered value, and compute its third moment. Orient by
the score's sign; retain both signs if its absolute value is at most `1e-12`.
Average predictions over the retained orientations (probabilities for binomial
outcomes), and tune one common penalty from their held-out losses. Recompute
the rule on the full training sample for refitting. Center relative to the
sample median before subtracting the residual mean to reduce cancellation.
This unsupervised rule uses one anchored fit on generic samples but can require
all anchors when the covariate distribution is empirically symmetric. Record
that cost and the selected orientation count. Evaluate it in a separate run
using the frozen primary pilot configuration so existing fits are reused.

Estimate the paired difference in integrated squared prediction error between
symmetric HAL and cross-validated selection of a single global anchor, averaged
over independent training samples from each specified data-generating process.
Use a common independent test sample and known conditional means. Report
absolute differences and risk ratios by scenario, with Monte Carlo uncertainty.

Compare anchored HAL, canonical orientation, selected anchor, uniformly
averaged anchor predictions, and symmetric HAL with the same zero-order basis
degree and knot counts. Use minimum-norm selection for every constituent fit.
For small dimension enumerate all sign patterns. Select the anchor and penalty
jointly by training-fold validation loss, with a declared deterministic tie rule.
Average all anchors with equal weights. Select one common penalty by the loss
of their averaged out-of-fold predictions; do not average separately tuned fits
or weight anchors using test outcomes.
An oracle chooses among the individually cross-validated anchor fits using the
known evaluation risk. It does not optimize jointly over anchors and penalties
and is unavailable to a practitioner. CART and random
forests provide secondary practical benchmarks.

Anchor selection itself need not be equivariant when distinct anchors tie.
Report that possibility rather than claiming invariance for a first-index tie
rule. Uniform averaging over all anchors is equivariant if each constituent
fit and the common tuning procedure transform consistently.

## First study: feasibility and mechanisms

Use two- and three-dimensional problems before increasing dimension. Include
an upper orthant, lower orthant, mixed orthant, interior box, two components
with conflicting preferred anchors, additive steps, and a smooth signal.
Add independent and correlated covariates, with a discrete/duplicate-covariate
stress case for solution selection. The conflicting-component example is the
main mechanism test; single-corner signals measure the price of flexibility.

Start with 20 training samples per cell, sample sizes 100 and 300, three knots
per coordinate, and signal-to-noise ratios 0.5 and 2. Use five common folds,
whole-training-sample outcome centering and standardization, and penalties
`0.003, 0.01, 0.03, 0.1, 0.3` on that scale. This scaling is fixed across
validation folds. Correlated covariates use a Gaussian copula with correlation
0.5. Discrete covariates are uniform on `{0, .25, .5, .75, 1}`, with the third
coordinate duplicating the second in dimension three.

Integrate squared error exactly for the uniform and discrete designs. For the
correlated design use a common independent sample of 10,000 covariate draws.
Also save error on that test sample for every design. CART and random forests
use the common test sample; their tuning follows the implementation defaults.
A mechanism comparison uses fixed knots `.25, .5, .75` for the conflicting
components at dimension two and SNR 2.
The 3,400 cases count scenario evaluations: fixed-knot cases reuse the matching
quantile-knot samples and folds. Signal variance is exact for independent and
discrete designs and estimated from a separate 200,000-point calibration sample
for the correlated design.
A nine-sample feasibility exercise preceded this pilot; its observed risks
are exploratory and do not constitute a test of comparative expected risk.
The present demonstrations and mathematical examples informed this design;
it is not a blinded discovery exercise.

For each fit save validation losses, selected penalty/anchor, held-out loss,
maximum prediction movement across codings, all optimality residuals, elapsed
method-specific fitting time including dictionary construction and tuning, and
failures. Record total elapsed time for each job separately: the shared suite
performs additional diagnostic fits that no single method requires. Store the
code revision, configuration, random seeds, and dependency versions with the
run. Never silently omit failed replicates. Report both the attempted sample
and any complete-pair analysis, including the reason each pair is missing.

Use paired means and Monte Carlo standard errors across independent training
samples. Display every planned scenario rather than pooling favorable ones.
Check whether any advantage survives equalizing tuning compute, and whether
it is driven by knot approximation or by optimization failures. For screening
larger grids, require zero uninvestigated fit failures and reflected prediction
movement below `1e-6 * max(1, sd(y))` on the designated test inputs. This is a
numerical acceptance criterion, not a theorem about arbitrary inputs.

## Decision gates

A method's central claim requires a proof under the declared transformation
family, a well-defined fitted function away from observed covariates, and
numerical checks of the complete tuned procedure and downstream inference.
No fixed percentage improvement in prediction risk is required. Report accuracy,
stability to perturbations, and computational cost as separate tradeoffs.

The completed prediction pilot gives mixed risk comparisons; retain these
findings without treating them as a failure of the invariance objective.
Discrete-case gains can reflect dictionary span. Before using such gains to
support a penalty mechanism, equalize the span with suitable support cutpoints.
Before expanding a population study, state its estimand, the mechanism-derived
conditions, a meaningful precision target, and the Monte Carlo uncertainty
needed to answer that question. Keep confirmatory draws independent of the
exploratory draws that informed the design.

## Theory and computation

The statistical collaborator can lead a literature and proof review of the
infimal-convolution formulation and strict gap relative to all single anchors.
Anchor dependence itself is established background; novelty must rest on a
precise result about mixed components, estimator definition, or computation.
Check whether the lattice entropy transfer adds enough beyond an application
of known convex-hull bounds to deserve a central claim.

In parallel with that scientific work, improve computation only after profiling
identifies the bottleneck. A sparse or matrix-free minimum-norm solver must
match primary objectives, full-dictionary secondary optimality, and off-sample
equivariance on the existing small problems. Training-column deduplication is
not an acceptable shortcut. Benchmark total time against enumerating anchors;
retain the dense implementation as a correctness benchmark during development,
without exposing competing public definitions of symmetric HAL.

## Causal extension

First isolate coding sensitivity on a fixed observational sample. Use the
Dehejia--Wahba treated sample and PSID comparison group, the same five outer
folds and ten inner folds, and three codings: original, continuous covariates
negated, and age in months with pretreatment earnings in thousands. Fit both
outcome regressions and the propensity score. Compare raw AIPW with AIPW and
ATT TMLE using common clipped outcomes and calibrated initial propensities.
Target on the pooled sample after cross-fitting the initial nuisances. Keep
the initial outcome clipping at `1e-6` separate from the numerical targeting
floor of `1e-12`, so initially clipped predictions can move along the path.
Check
the empirical influence-function mean and stability across targeting step
sizes; do not equate a software convergence flag with a solved score equation.

Separate package-default behavior from solver and scaling effects. Repeat
with a tighter convergence threshold and a larger iteration budget, with
training-sample covariate standardization, and with both controls. Preserve
all cases and record warnings by fold and nuisance component. A change under
recoding establishes sensitivity, not which estimate is correct.

The completed controls retain the point-estimate sign reversal under each
setting. Standardization removes observed unit dependence to numerical
precision with either solver setting; tighter optimization alone does not.
Convergence warnings persist, so these controls do not establish that every
nuisance validation path completed or that optimization effects are absent.

Apply the canonical remedy directly to the first-order package as well. For
each outer-training nuisance sample, standardize covariates using its own
means and standard deviations, then apply the same normalized third-moment
orientation rule, averaging all ambiguous directions. Each orientation uses
the original package's basis defaults and penalty selection with shared inner
folds. Average its separately tuned response predictions. This differs from
the controlled zero-order comparator, which rebuilds orientation inside every
inner training fold and tunes one common penalty for the averaged predictor.
The first-order preprocessing is learned from outer-training covariates only;
held-out outer-fold covariates and outcomes are not used to choose directions.
Compare this procedure with the existing scaled first-order control to assess
what adding orientation contributes. Check nuisance predictions and final ATT
estimates under the original, reflected, and different-units codings.

Then use a known-effect semi-synthetic design, common outer and inner folds,
and identical propensity clipping across methods. Report bias, RMSE, interval
coverage with Monte Carlo uncertainty, nuisance loss, overlap diagnostics, and
coding movement. Separate orientation, knot placement, and scaling. The
LaLonde data can illustrate an application but proximity to its experimental
estimate cannot validate the estimator or replace the known-effect study.

The deliverable is a reproducible account of causal coding sensitivity and
comparisons of practical remedies, including their numerical reliability,
prediction risk, causal accuracy, and computational cost. The public paper
should contain only claims supported by the mathematics and the completed,
reproducible experiments.


### Known-effect design

Freeze the empirical distribution of the 2,675 treated/PSID covariate rows.
Fit a logistic treatment model with age and its square, education, demographic
indicators, pretreatment earnings, and zero-earnings indicators; cap the true
propensity at 0.95. Fit the untreated mean to PSID observations using natural
splines in age and pretreatment earnings plus the other covariates. Save the
resulting true probabilities and means, regression coefficients, residual
standard deviation, and source row identifiers.

For each replicate resample 2,675 covariate rows, draw treatment from the true
propensity, and draw a Gaussian outcome with untreated mean plus 1,800 times
treatment. These are earnings-like continuous outcomes and can be negative.
The ATT is exactly 1,800; the simulation population is conditional on the fixed
source dataset. Compute the nonparametric ATT efficiency benchmark once from
that entire population, not from each resampled design. Use observed-range
outcome scaling for TMLE; do not truncate simulated outcomes to be nonnegative.

Run two draws for computational feasibility, then 20 for an exploratory pilot.
The pilot uses draws 1 through 20, including the two feasibility draws; it is
not an independent validation sample.
Compare true initial nuisances, anchored zero-order HAL, canonical orientation,
symmetric zero-order HAL, first-order hal9001 defaults, and the canonicalized
first-order wrapper, under original and continuous-negated
codings. Reuse each draw and every fold assignment. Gaussian zero-order nuisance
fits standardize outcomes within the outer training sample, use penalties
`.003, .01, .03, .1`, and restore prediction units. Binomial penalties are
`.002, .01, .05, .2`. Cross-validation loss ties use the normalized outcome scale.
The first-order package benchmark follows its own defaults.

Report raw oracle AIPW separately from oracle-initialized TMLE: the latter's
clipping and pooled calibration alter its true initial nuisance inputs. Preserve
both outcome regressions and the propensity predictions. Assess untreated-mean
risk both over all covariates and weighted by true propensity for the treated
population, alongside propensity risk and a plug-in nuisance-error drift
calculation. Report bias, RMSE, paired coding movement, Monte Carlo standard
errors, and exact binomial intervals for empirical coverage. Twenty replicates
cannot establish calibrated coverage. Use pilot uncertainty to size any further
accuracy study before drawing additional replicates.

Compare canonical first-order HAL with its package-default counterpart,
canonical zero-order HAL with anchored HAL, and symmetric HAL with both
anchored and canonical zero-order HAL. Pair by replicate, coding, and causal
estimator. For errors `e` and `r`, report the mean of `e^2 - r^2` with MCSE
`sd(e^2 - r^2) / sqrt(R)`. For the ratio of RMSEs, let `m = mean(e^2)` and
`v = mean(r^2)`. When both are positive, the delta-method MCSE on the log scale
is `sd(e^2/m - r^2/v) / (2*sqrt(R))`. Use a pointwise approximate interval with
the `t` critical value on `R - 1` degrees of freedom, exponentiating the log
ratio interval. These are exploratory Monte Carlo intervals, not confidence
intervals for the causal effect or multiplicity-adjusted claims. When empirical
MSE is zero, leave log-ratio uncertainty undefined rather than divide by zero.

Require the complete replicate-by-method-by-coding-by-estimator grid and the
same frozen population model when collecting multiple run directories. Report
absolute standard-error movement alongside ATT movement. Count clipping and
nuisance warnings once per fitted replicate and coding. Preparation clipping
describes inputs shared by calibrated AIPW and TMLE, not raw AIPW; nuisance
warnings do not include targeting messages. Wall times under concurrent load
are descriptive and cannot establish an isolated computational advantage.

After collecting all six methods for pilot draws 1 through 20, run
`make semisynth-paper`. This recomputes the statistical summaries and generates
separate accuracy tables for raw AIPW, prepared AIPW, and TMLE, a paired coding
movement table, and the planned paired RMSE comparisons. Accuracy tables show
both codings, Monte Carlo uncertainty, exact coverage intervals, and counts of
draws with nuisance warnings. Movement remains in dollars with scientific
notation for small nonzero values; accuracy is in thousands of dollars.
The paper-table target rejects an incomplete method or replicate set. The
completed 20-draw pilot provides the manuscript comparisons; the two-draw
feasibility subset is not a separate validation study. Canonical procedures
agree exactly across sign codings and symmetric HAL agrees to numerical
precision. Accuracy rankings depend on the causal estimator. These findings
support the invariance objective without establishing general risk or coverage
advantages. Any additional population study still requires a mechanism-derived
question and precision calculation before new draws.
