# ==============================================================================
# Generate Legacy Fit Objects for Backward Compatibility Testing
# ==============================================================================
#
# This script installs old versions of bgms from Posit Package Manager binary
# snapshots and creates minimal fit objects for testing backward compatibility
# of extractor functions.
#
# Run this script manually when needed to regenerate fixtures.
#
# Output: tests/testthat/fixtures/legacy/
#   - fit_v{VERSION}.rds for each historical version
#
# These fixtures are NOT shipped with the package (tests/ is excluded from
# installed packages). Legacy tests skip on CRAN where fixtures aren't available.
#
# ==============================================================================

library(callr)

# Output directory (in tests, not inst)
legacy_dir <- file.path("tests", "testthat", "fixtures", "legacy")
dir.create(legacy_dir, recursive = TRUE, showWarnings = FALSE)

# Minimal data for fitting
set.seed(42)
test_data <- matrix(sample(0:1, 60, replace = TRUE), ncol = 3)
colnames(test_data) <- c("X1", "X2", "X3")

# Version → snapshot date mapping (use date when version was current)
# Dates derived from CRAN archive release dates
version_snapshots <- list(
  "0.1.3"   = "2024-03-15",   # Released 2024-02-25
  "0.1.3.1" = "2024-06-01",   # Released 2024-05-15
  "0.1.4"   = "2024-11-01",   # Released 2024-10-21
  "0.1.4.1" = "2024-11-20",   # Released 2024-11-12
  "0.1.4.2" = "2024-12-15",   # Released 2024-12-05
  "0.1.6.0" = "source",       # Install from source (no binary available)
  "0.1.6.1" = "source"        # Install from source (no binary available)
)

# Helper function to install and run in isolated environment
create_legacy_fit <- function(version, snapshot_date, output_name) {
  cat(sprintf("\n=== Creating fixture: %s (v%s) ===\n", output_name, version))
  
  tmp_lib <- tempfile("bgms_lib_")
  dir.create(tmp_lib)
  
  tryCatch({
    if (snapshot_date == "source") {
      # Install from CRAN archive source
      url <- paste0("https://cran.r-project.org/src/contrib/Archive/bgms/bgms_", version, ".tar.gz")
      cat("Installing bgms", version, "from source:", url, "\n")
      install.packages(
        url,
        repos = NULL,
        type = "source",
        lib = tmp_lib,
        quiet = TRUE
      )
    } else {
      # Install from Posit Package Manager binary snapshot
      repos <- paste0("https://packagemanager.posit.co/cran/", snapshot_date)
      cat("Installing bgms", version, "from", repos, "(binary)\n")
      install.packages(
        "bgms",
        lib = tmp_lib,
        repos = repos,
        type = "binary",
        quiet = TRUE
      )
    }
    
    # Verify installation
    if (!"bgms" %in% list.files(tmp_lib)) {
      stop("bgms was not installed")
    }
    
    # Run fit in separate R process with isolated library
    result <- callr::r(
      function(data, lib_path) {
        .libPaths(c(lib_path, .libPaths()))
        library(bgms, lib.loc = lib_path)
        pkg_version <- as.character(packageVersion("bgms"))
        cat("bgms version:", pkg_version, "\n")
        
        # Create fit with edge selection (save=TRUE for raw samples)
        fit <- bgm(data, iter = 100, burnin = 50, save = TRUE)
        
        # Return fit object
        fit
      },
      args = list(data = test_data, lib_path = tmp_lib),
      show = TRUE
    )
    
    # Save fixture
    output_path <- file.path(legacy_dir, paste0(output_name, ".rds"))
    saveRDS(result, output_path)
    cat("Saved:", output_path, "\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    return(FALSE)
  }, finally = {
    unlink(tmp_lib, recursive = TRUE)
  })
}

# ==============================================================================
# Generate fixtures for all versions
# ==============================================================================

results <- list()

for (version in names(version_snapshots)) {
  snapshot_date <- version_snapshots[[version]]
  output_name <- paste0("fit_v", version)
  
  success <- create_legacy_fit(
    version = version,
    snapshot_date = snapshot_date,
    output_name = output_name
  )
  
  results[[version]] <- success
}

# Summary
cat("\n=== Summary (bgm) ===\n")
for (version in names(results)) {

  status <- if (results[[version]]) "SUCCESS" else "FAILED"
  cat(sprintf("  %s: %s\n", version, status))
}

# ==============================================================================
# Generate bgmCompare fixtures
# ==============================================================================
#
# bgmCompare was introduced in version 0.1.4. Generate fixtures for versions
# 0.1.4+ to test backward compatibility of bgmCompare extractor functions.
#
# ==============================================================================

# Create two-group test data for bgmCompare
set.seed(42)
group1_data <- matrix(sample(0:1, 60, replace = TRUE), ncol = 3)
group2_data <- matrix(sample(0:1, 60, replace = TRUE), ncol = 3)
colnames(group1_data) <- c("X1", "X2", "X3")
colnames(group2_data) <- c("X1", "X2", "X3")

# Versions that have bgmCompare (introduced in 0.1.4)
bgmcompare_versions <- list(
  "0.1.4"   = "2024-11-01",
  "0.1.4.1" = "2024-11-20",
  "0.1.4.2" = "2024-12-15",
  "0.1.6.0" = "source",
  "0.1.6.1" = "source"
)

# Helper function to create bgmCompare legacy fit
create_legacy_bgmcompare_fit <- function(version, snapshot_date, output_name) {
  cat(sprintf("\n=== Creating bgmCompare fixture: %s (v%s) ===\n", output_name, version))
  
  tmp_lib <- tempfile("bgms_lib_")
  dir.create(tmp_lib)
  
  tryCatch({
    if (snapshot_date == "source") {
      # Install from CRAN archive source
      url <- paste0("https://cran.r-project.org/src/contrib/Archive/bgms/bgms_", version, ".tar.gz")
      cat("Installing bgms", version, "from source:", url, "\n")
      install.packages(
        url,
        repos = NULL,
        type = "source",
        lib = tmp_lib,
        quiet = TRUE
      )
    } else {
      # Install from Posit Package Manager binary snapshot
      repos <- paste0("https://packagemanager.posit.co/cran/", snapshot_date)
      cat("Installing bgms", version, "from", repos, "(binary)\n")
      install.packages(
        "bgms",
        lib = tmp_lib,
        repos = repos,
        type = "binary",
        quiet = TRUE
      )
    }
    
    # Verify installation
    if (!"bgms" %in% list.files(tmp_lib)) {
      stop("bgms was not installed")
    }
    
    # Run fit in separate R process with isolated library
    result <- callr::r(
      function(data_x, data_y, lib_path) {
        .libPaths(c(lib_path, .libPaths()))
        library(bgms, lib.loc = lib_path)
        pkg_version <- as.character(packageVersion("bgms"))
        cat("bgms version:", pkg_version, "\n")
        
        # Create bgmCompare fit with difference_selection (save=TRUE for raw samples)
        # Old API uses separate x and y data frames
        fit <- bgmCompare(
          x = data_x,
          y = data_y,
          iter = 100,
          burnin = 50,
          difference_selection = TRUE,
          save = TRUE
        )
        
        # Return fit object
        fit
      },
      args = list(data_x = group1_data, data_y = group2_data, lib_path = tmp_lib),
      show = TRUE
    )
    
    # Save fixture
    output_path <- file.path(legacy_dir, paste0(output_name, ".rds"))
    saveRDS(result, output_path)
    cat("Saved:", output_path, "\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    return(FALSE)
  }, finally = {
    unlink(tmp_lib, recursive = TRUE)
  })
}

# Generate bgmCompare fixtures
bgmcompare_results <- list()

for (version in names(bgmcompare_versions)) {
  snapshot_date <- bgmcompare_versions[[version]]
  output_name <- paste0("bgmcompare_v", version)
  
  success <- create_legacy_bgmcompare_fit(
    version = version,
    snapshot_date = snapshot_date,
    output_name = output_name
  )
  
  bgmcompare_results[[version]] <- success
}

# Summary
cat("\n=== Summary (bgmCompare) ===\n")
for (version in names(bgmcompare_results)) {
  status <- if (bgmcompare_results[[version]]) "SUCCESS" else "FAILED"
  cat(sprintf("  %s: %s\n", version, status))
}

cat("\nAll fixtures created in:", legacy_dir, "\n")
