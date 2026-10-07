# Installs the R package nebula and the dependency versions it needs.
# Run with:  Rscript install_nebula.R
repos <- "https://cloud.r-project.org"
# nebula's declared minimum versions are too low for the API it uses; upgrade these if older
need <- function(p, v) !requireNamespace(p, quietly = TRUE) || packageVersion(p) < v
up <- c(parallelly = "1.38.0", future = "1.40.0", future.apply = "1.11.0", globals = "0.16.3", doFuture = "1.0.0")
todo <- names(up)[mapply(need, names(up), up)]
if (length(todo)) install.packages(todo, repos = repos)
if (!requireNamespace("nebula", quietly = TRUE)) {
  ok <- tryCatch({install.packages("nebula", repos = repos); TRUE}, error = function(e) FALSE)
  if (!requireNamespace("nebula", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
    remotes::install_github("lhe17/nebula", upgrade = "never")
  }
}
ok <- tryCatch({
  suppressPackageStartupMessages(library(nebula)); data(sample_data, package = "nebula")
  re <- nebula(sample_data$count, sample_data$sid, pred = model.matrix(~X1, data = sample_data$pred),
               ncore = 1, verbose = FALSE)
  nrow(re$summary) > 0
}, error = function(e) {message("smoke test failed: ", conditionMessage(e)); FALSE})
cat("NEBULA_OK", ok, "\n")
