# ===========================================================================
# Mixed MRF Validation Suite — Master Runner
# ===========================================================================
# Runs all 7 test groups sequentially. Each group produces:
#   - Console output with numerical summaries
#   - PDF plots in dev/tests/validation/output/
#
# Usage:
#   cd /path/to/bgms
#   Rscript dev/tests/validation/run_all.R
#
# Or run individual groups:
#   source("dev/tests/validation/group1_parameter_recovery.R")
#
# Prerequisites:
#   - bgms installed or load_all'able
#   - mixedGM installed (for group 4)
# ===========================================================================

t0 = proc.time()
cat("===================================================================\n")
cat("  Mixed MRF Validation Suite\n")
cat(sprintf("  Started: %s\n", Sys.time()))
cat("===================================================================\n\n")

# Ensure output directory exists
out_dir = file.path("dev", "tests", "validation", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Group 1: Parameter Recovery
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 1: Parameter Recovery                                     #\n")
cat("###################################################################\n")
tryCatch(
  source("dev/tests/validation/group1_parameter_recovery.R", local = TRUE),
  error = function(e) cat("GROUP 1 FAILED:", conditionMessage(e), "\n")
)

# Group 2: Metropolis vs NUTS
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 2: Metropolis vs NUTS                                     #\n")
cat("###################################################################\n")
tryCatch(
  source("dev/tests/validation/group2_mh_vs_nuts.R", local = TRUE),
  error = function(e) cat("GROUP 2 FAILED:", conditionMessage(e), "\n")
)

# Group 3: Conditional vs Marginal PL
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 3: Conditional vs Marginal PL                             #\n")
cat("###################################################################\n")
tryCatch(
  source("dev/tests/validation/group3_cond_vs_marg.R", local = TRUE),
  error = function(e) cat("GROUP 3 FAILED:", conditionMessage(e), "\n")
)

# Group 4: Cross-package (bgms vs mixedGM)
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 4: bgms vs mixedGM Cross-Validation                      #\n")
cat("###################################################################\n")
if(requireNamespace("mixedGM", quietly = TRUE)) {
  tryCatch(
    source("dev/tests/validation/group4_cross_package.R", local = TRUE),
    error = function(e) cat("GROUP 4 FAILED:", conditionMessage(e), "\n")
  )
} else {
  cat("SKIPPED: mixedGM package not available.\n")
}

# Group 5: Edge Detection Accuracy
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 5: Edge Detection Accuracy                                #\n")
cat("###################################################################\n")
tryCatch(
  source("dev/tests/validation/group5_edge_detection.R", local = TRUE),
  error = function(e) cat("GROUP 5 FAILED:", conditionMessage(e), "\n")
)

# Group 6: MCMC Diagnostics
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 6: MCMC Diagnostics                                       #\n")
cat("###################################################################\n")
tryCatch(
  source("dev/tests/validation/group6_diagnostics.R", local = TRUE),
  error = function(e) cat("GROUP 6 FAILED:", conditionMessage(e), "\n")
)

# Group 7: Empirical Data Smoke Tests
cat("\n\n")
cat("###################################################################\n")
cat("# GROUP 7: Empirical Data Smoke Tests                             #\n")
cat("###################################################################\n")
tryCatch(
  source("dev/tests/validation/group7_empirical.R", local = TRUE),
  error = function(e) cat("GROUP 7 FAILED:", conditionMessage(e), "\n")
)

# ===========================================================================
# Final summary
# ===========================================================================
elapsed = (proc.time() - t0)["elapsed"]
cat("\n\n")
cat("===================================================================\n")
cat(sprintf("  Validation suite complete. Total time: %.1f minutes\n",
            elapsed / 60))
cat(sprintf("  Finished: %s\n", Sys.time()))
cat(sprintf("  Output directory: %s\n", normalizePath(out_dir)))
cat("===================================================================\n")

# List generated PDFs
pdfs = list.files(out_dir, pattern = "\\.pdf$", full.names = FALSE)
if(length(pdfs) > 0) {
  cat("\n  Generated plots:\n")
  for(f in pdfs) cat(sprintf("    %s\n", f))
}
