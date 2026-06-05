# Quiet bgms's advisory verbose output (warmup-length warnings, data-cleaning
# messages) during the test run. Tests that specifically assert this output set
# options(bgms.verbose = TRUE) locally, which overrides this default.
options(bgms.verbose = FALSE)
