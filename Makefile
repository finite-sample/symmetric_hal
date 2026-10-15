.PHONY: test lint check sims recode report paper clean

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

# Figures and tables from data/ alone.
report:
	Rscript R/fig_coding.R
	Rscript R/sim_report.R tier1
	Rscript R/recode_report.R

paper: report
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error ms.tex

clean:
	cd paper && latexmk -C ms.tex
