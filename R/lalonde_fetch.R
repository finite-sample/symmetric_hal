# Rscript R/lalonde_fetch.R
# Writes data/lalonde_dw.csv: the Dehejia-Wahba NSW subsample (185 treated,
# 260 experimental controls) and the PSID-1 comparison group (2,490), from
# Dehejia's NBER page. `group` is "treated", "experimental" or "psid".
suppressMessages(library(tidyverse))
base <- "https://users.nber.org/~rdehejia/data/"
vars <- c("treat", "age", "educ", "black", "hisp", "married", "nodegree",
          "re74", "re75", "re78")
files <- c(treated = "nswre74_treated.txt",
           experimental = "nswre74_control.txt", psid = "psid_controls.txt")
imap(files, \(f, group) {
  read_table(paste0(base, f), col_names = vars, show_col_types = FALSE) |>
    mutate(group = group, .before = 1)
}) |>
  list_rbind() |>
  write_csv("data/lalonde_dw.csv")
