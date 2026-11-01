export R_LIBS := $(CURDIR)/.Rlib:$(R_LIBS)

.PHONY: test lint check example sims recode lalonde semisynth audit data report benchmark-reports paper clean pilot causal causal-report causal-sensitivity
.PHONY: oracle-report pilot-report causal-control-report semisynth-paper

CORES ?= 8
REPS ?= 50
PILOT_DIR ?= results_current_pilot
CANONICAL_DIR ?= results_current_canonical
CAUSAL_DIR ?= results_current_causal
SEMISYNTH_REPS ?= 20

test:
	Rscript -e 'testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'

lint:
	Rscript -e 'l <- lintr::lint_dir("."); print(l); quit(status = as.integer(length(l) > 0))'

check: lint test

example:
	Rscript R/example.R
	Rscript R/fig_coding.R

pilot:
	mkdir -p $(PILOT_DIR)
	Rscript R/pilot_run.R $(PILOT_DIR) pilot $(CORES) 20 >> $(PILOT_DIR)/run.log 2>&1
	mkdir -p $(CANONICAL_DIR)
	Rscript R/canonical_run.R $(PILOT_DIR) $(CANONICAL_DIR) $(CORES) >> $(CANONICAL_DIR)/run.log 2>&1
	Rscript R/pilot_collect.R $(PILOT_DIR) $(CANONICAL_DIR)
	Rscript R/pilot_report.R

causal:
	mkdir -p $(CAUSAL_DIR)
	Rscript R/causal_run.R $(CAUSAL_DIR) $(CORES) >> $(CAUSAL_DIR)/run.log 2>&1
	Rscript R/causal_collect.R $(CAUSAL_DIR)
	Rscript R/causal_report.R

causal-report:
	Rscript R/causal_report.R

causal-control-report:
	Rscript R/causal_control_report.R

causal-sensitivity:
	Rscript R/causal_sensitivity.R

# Long runs. Each fit is checkpointed in results_current_*/ (gitignored), so a stopped
# run resumes; collecting writes the small CSVs in data/ that are committed.
sims:
	mkdir -p results_current_tier1
	Rscript R/sim_run.R tier1 $(REPS) results_current_tier1 $(CORES) >> results_current_tier1/run.log 2>&1
	Rscript R/sim_collect.R results_current_tier1 tier1

recode:
	mkdir -p results_current_recode
	Rscript R/recode_run.R results_current_recode 21 $(CORES) >> results_current_recode/run.log 2>&1
	Rscript R/recode_collect.R results_current_recode

lalonde:
	mkdir -p results_current_lalonde
	Rscript R/lalonde_run.R results_current_lalonde 21 $(CORES) >> results_current_lalonde/run.log 2>&1
	Rscript R/lalonde_collect.R results_current_lalonde
	mkdir -p results_current_lalonde_units
	Rscript R/lalonde_units_run.R results_current_lalonde_units $(CORES) >> results_current_lalonde_units/run.log 2>&1
	Rscript R/lalonde_collect.R results_current_lalonde_units data/lalonde_units.csv

semisynth:
	mkdir -p results_current_semisynth
	Rscript R/semisynth_run.R results_current_semisynth $(SEMISYNTH_REPS) $(CORES) >> results_current_semisynth/run.log 2>&1
	SEMISYNTH_CORES=$(CORES) Rscript R/semisynth_collect.R results_current_semisynth
	Rscript R/semisynth_report.R

audit:
	mkdir -p results_current_audit
	Rscript R/audit_run.R results_current_audit $(CORES) >> results_current_audit/run.log 2>&1
	Rscript R/audit_collect.R results_current_audit

# Re-download the Dehejia-Wahba and PSID-1 files into data/lalonde_dw.csv.
data:
	Rscript R/lalonde_fetch.R

oracle-report:
	Rscript R/semisynth_oracle_report.R

pilot-report:
	Rscript R/pilot_report.R

semisynth-paper:
	Rscript R/semisynth_paper.R

report: example causal-report causal-control-report causal-sensitivity oracle-report pilot-report semisynth-paper

benchmark-reports:
	Rscript R/sim_report.R tier1
	Rscript R/recode_report.R
	Rscript R/lalonde_report.R
	Rscript R/semisynth_report.R
	Rscript R/audit_report.R

paper: report
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error ms.tex

clean:
	cd paper && latexmk -C ms.tex
