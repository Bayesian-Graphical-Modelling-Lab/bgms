Rcpp::sourceCpp("src/progress_manager.cpp")
runMCMC_parallel(3, 500)
runMCMC_parallel(3, 500, useUnicode = FALSE)
runMCMC_parallel(3, 500, useUnicode = FALSE, display_progress = TRUE)