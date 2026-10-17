.PHONY: test lint check sims recode lalonde semisynth audit data report paper clean

CORES ?= 8
REPS ?= 50

test:
	Rscript -e 'testthat::test_dir("tests/testthat", stop_on_failure = TRUE)'

lint:
	Rscript -e 'l <- lintr::lint_dir("."); print(l); quit(status = as.integer(length(l) > 0))'

check: lint test

# Long runs. Each fit is checkpointed in results_*/ (gitignored), so a stopped
# run resumes; collecting writes the small CSVs in data/ that are committed.
sims:
	mkdir -p results_tier1
	Rscript R/sim_run.R tier1 $(REPS) results_tier1 $(CORES) >> results_tier1/run.log 2>&1
	Rscript R/sim_collect.R results_tier1 tier1

recode:
	mkdir -p results_recode
	Rscript R/recode_run.R results_recode 21 $(CORES) >> results_recode/run.log 2>&1
	Rscript R/recode_collect.R results_recode

lalonde:
	mkdir -p results_lalonde
	Rscript R/lalonde_run.R results_lalonde 21 $(CORES) >> results_lalonde/run.log 2>&1
	Rscript R/lalonde_collect.R results_lalonde
	mkdir -p results_lalonde_units
	Rscript R/lalonde_units_run.R results_lalonde_units $(CORES) >> results_lalonde_units/run.log 2>&1
	Rscript R/lalonde_collect.R results_lalonde_units data/lalonde_units.csv

semisynth:
	mkdir -p results_semisynth
	Rscript R/semisynth_run.R results_semisynth 66 $(CORES) >> results_semisynth/run.log 2>&1
	Rscript R/semisynth_collect.R results_semisynth

audit:
	mkdir -p results_audit
	Rscript R/audit_run.R results_audit $(CORES) >> results_audit/run.log 2>&1
	Rscript R/audit_collect.R results_audit

# Re-download the Dehejia-Wahba and PSID-1 files into data/lalonde_dw.csv.
data:
	Rscript R/lalonde_fetch.R

# Figures and tables from data/ alone.
report:
	Rscript R/fig_coding.R
	Rscript R/sim_report.R tier1
	Rscript R/recode_report.R
	Rscript R/lalonde_report.R
	Rscript R/semisynth_report.R
	Rscript R/audit_report.R

paper: report
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error ms.tex

clean:
	cd paper && latexmk -C ms.tex
