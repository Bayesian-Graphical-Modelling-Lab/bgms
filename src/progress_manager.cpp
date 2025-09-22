#include "progress_manager.h"

ProgressManager::ProgressManager(int nChains_, int nIter_, int nWarmup_, int printEvery_, int progress_type_, bool useUnicode_)
  : nChains(nChains_), nIter(nIter_ + nWarmup_), nWarmup(nWarmup_), progress(nChains_), printEvery(printEvery_),
    progress_type(progress_type_), useUnicode(useUnicode_) { // +2 for total and time lines

  for (int i = 0; i < nChains; i++) progress[i] = 0;
  start = Clock::now();
  lastPrint = Clock::now();

  // Check if we're in RStudio
  Rcpp::Environment base("package:base");
  Rcpp::Function getOption = base["getOption"];
  Rcpp::Function Sysgetenv = base["Sys.getenv"];
  SEXP s = Sysgetenv("RSTUDIO");
  isRStudio = Rcpp::as<std::string>(s) == "1";

  no_spaces_for_total = 3 + static_cast<int>(std::log10(nChains));
  if (progress_type == 1) no_spaces_for_total = 0; // no total line
  total_padding = std::string(no_spaces_for_total, ' ');

  if (isRStudio) {
    consoleWidth = getConsoleWidth();
    lineWidth = std::max(10, std::min(consoleWidth - 25, 70));
  } else {
    // For terminal, use default console width
    consoleWidth = 80;
    lineWidth = 70;
  }

  // cleverly determine barwidth so no line wrapping occurs
  if (lineWidth <= 5) {
    // TODO: we don't want to print anything in this case
    barWidth = 0;
  } else if (lineWidth < 20) {
    barWidth = lineWidth - 10;
  } else if (lineWidth < 40) {
    barWidth = lineWidth - 15;
  } else { // > 40
    barWidth = std::min(40, lineWidth - 30);
  }

  if (isRStudio) {
    barWidth = std::max(10, barWidth - 20); // minimum bar width of 10 for RStudio
  }

  // Set up theme
  setupTheme();
  update_prefixes(consoleWidth);
}

void ProgressManager::update(int chainId) {
  progress[chainId]++;

  // Only chain 0 actually does the printing/ checking for user interrupts
  if (chainId != 0) return;

  if (progress[chainId] % printEvery == 0) {
    auto now = Clock::now();
    std::chrono::duration<double> sinceLast = now - lastPrint;

    // Throttle printing to avoid spamming
    if (progress_type != 0 && sinceLast.count() >= 0.5) {
      print();
      lastPrint = now;
    }
  }

  // Check for user interrupts and console width changes less frequently to reduce overhead
  if (chainId == 0 && progress[chainId] % (printEvery * 5) == 0) {
    needsToExit = checkInterrupt();
    // Also check for console width changes occasionally
    checkConsoleWidthChange();
  }
}

void ProgressManager::finish() {

  if (progress_type == 0) return; // No progress display

  // Mark all chains as complete and print one final time
  for (int i = 0; i < nChains; i++) {
    progress[i] = nIter;
  }
  print();

}

bool ProgressManager::shouldExit() const {
  return needsToExit;
}

void ProgressManager::checkConsoleWidthChange() {
  if (!isRStudio) return;

  Rcpp::Environment base("package:base");
  Rcpp::Function getOption = base["getOption"];
  int currentWidth = getConsoleWidth();

  if (prevConsoleWidth == -1) {
    // First time, just store the current width
    prevConsoleWidth = consoleWidth;
    widthChanged = false;
    return;
  }

  if (currentWidth != consoleWidth && currentWidth > 0) {
    // Width has changed
    prevConsoleWidth = consoleWidth;
    consoleWidth = currentWidth;
    widthChanged = true;
  } else {
    widthChanged = false;
  }
}

int ProgressManager::getConsoleWidth() {
  Rcpp::Environment base("package:base");
  Rcpp::Function getOption = base["getOption"];
  SEXP s = getOption("width", 0);
  int width = Rcpp::as<int>(s);
  // Remove the +3 adjustment to test actual console width
  return width + 3;
}

std::string ProgressManager::formatProgressBar(int chainId, int current, int total, double fraction, bool isTotal) {
    std::ostringstream builder;

    double exactFilled = fraction * barWidth;
    int filled = std::max(0, std::min(int(exactFilled), barWidth));

    // Build progress bar with theme
    std::string progressBar = lhsToken;

    // Add filled tokens
    for (int i = 0; i < filled; i++) {
      progressBar += filledToken;
    }

    // Add partial token if needed
    if (filled < barWidth) {
      double partialAmount = exactFilled - filled;
      if (partialAmount > 0) {
        if (partialAmount > 0.5) {
          progressBar += partialTokenMore;
        } else {
          progressBar += partialTokenLess;
        }
        filled++; // Account for the partial token
      }
    }

    // Add empty tokens
    for (int i = filled; i < barWidth; i++) {
      progressBar += emptyToken;
    }

    progressBar += rhsToken;

    if (isTotal) {

        std::string warmupOrSampling = isWarmupPhase() ? "(Warmup)" : "(Sampling)";
        builder << total_prefix << total_padding << warmupOrSampling << ": " << progressBar << " " << current << "/" << total
                << " (" << std::fixed << std::setprecision(1) << fraction * 100 << "%)";

    } else {

        std::string warmupOrSampling = isWarmupPhase(chainId - 1) ? " (Warmup)" : " (Sampling)";
        builder << chain_prefix << " " << chainId << warmupOrSampling << ": " << progressBar << " " << current << "/" << total
                << " (" << std::fixed << std::setprecision(1) << fraction * 100 << "%)";
    }

    return builder.str();
}

std::string ProgressManager::formatTimeInfo(int elapsed, int eta) {
  std::ostringstream builder;
  builder << "Elapsed: " << elapsed << "s | ETA: " << eta << "s";
  return builder.str();
}

void ProgressManager::setupTheme() {
  if (useUnicode) {
    // Unicode theme
    lhsToken = "❨";
    rhsToken = "❩";
    // filledToken = "\x1b[34m━\x1b[0m";  // Blue filled
    // emptyToken = "\x1b[37m━\x1b[0m";   // Gray empty
    // partialTokenMore = "\x1b[34m╸\x1b[0m";  // Blue partial (> 0.5)
    // partialTokenLess = "\x1b[37m╺\x1b[0m";  // Gray partial (< 0.5)
    filledToken = "━";  // Blue filled
    emptyToken = " ";   // Gray empty
    partialTokenMore = "╸";  // Blue partial (> 0.5)
    partialTokenLess = " ";  // Gray partial (< 0.5)
  } else {
    // Classic theme
    lhsToken = "[";
    rhsToken = "]";
    filledToken = "=";
    emptyToken = " ";
    partialTokenMore = " ";
    partialTokenLess = " ";
  }
}

size_t ProgressManager::getVisualLength(const std::string& str) {
  size_t visualLength = 0;
  bool inEscapeSequence = false;

  for (size_t i = 0; i < str.length(); i++) {
    if (str[i] == '\x1b' && i + 1 < str.length() && str[i + 1] == '[') {
      inEscapeSequence = true;
      i++; // Skip the '['
    } else if (inEscapeSequence && str[i] == 'm') {
      inEscapeSequence = false;
    } else if (!inEscapeSequence) {
      visualLength++;
    }
  }

  return visualLength;
}

void ProgressManager::print() {
  std::lock_guard<std::mutex> lock(printMutex);

  auto now = Clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();

  int totalWork = nChains * nIter;
  int done = std::reduce(progress.begin(), progress.end());
  double fracTotal = double(done) / totalWork;
  // should actually be the eta of the slowest chain!
  double eta = (fracTotal > 0) ? elapsed / fracTotal - elapsed : 0.0;

  std::ostringstream out;
  int totalChars = 0;
  int lineIndex = 0;

  // if this is not the first print, delete previous content
  if (lastPrintedChars > 0) {

    if (isRStudio) {
      out << "\x1b[" << std::to_string(lastPrintedChars) << "D";
    } else {
      // Move cursor up to start of our content and clear everything
      for (int i = 0; i < lastPrintedLines; i++) {
        out << "\x1b[1A\x1b[2K"; // Move up one line and clear entire line
      }
    }
  }

  if (progress_type == 2) {
    // Print progress for each chain
    for (int i = 0; i < nChains; i++) {
      double frac = double(progress[i]) / nIter;
      std::string chainProgress = formatProgressBar(i + 1, progress[i], nIter, frac);
      maybePadToLength(chainProgress);
      out << chainProgress << "\n";
      totalChars += chainProgress.length() + 1; // +1 for newline
    }

    // Print total progress
    std::string totalProgress = formatProgressBar(0, done, totalWork, fracTotal, true);
    maybePadToLength(totalProgress);
    out << totalProgress << "\n";
    totalChars += totalProgress.length() + 1; // +1 for newline

    // Print time info
    std::string timeInfo = formatTimeInfo(int(elapsed), int(eta));
    maybePadToLength(timeInfo);
    out << timeInfo << "\n";
    totalChars += timeInfo.length() + 1; // +1 for newline

    // Track total lines printed (chains + total + time)
    lastPrintedLines = nChains + 2; // used in a generic terminal
    lastPrintedChars = totalChars;  // used by RStudio

  } else if (progress_type == 1) {

    // Print total progress
    std::string totalProgress = formatProgressBar(0, done, totalWork, fracTotal, true);
    maybePadToLength(totalProgress);

    // Print time info
    totalProgress += " " + formatTimeInfo(int(elapsed), int(eta));
    maybePadToLength(totalProgress);

    if (done < totalWork) {
      out << totalProgress << "\r";
    } else {
      out << totalProgress << "\n";
    }

    // we do not set lastPrintedChars or lastPrintedLines here since we always overwrite the same line

  }

  Rcpp::Rcout << out.str();

}

void ProgressManager::update_prefixes(int width) {
  if (width < 20) {
    chain_prefix = "C";
    total_prefix = "T";
  } else if (width < 30) {
    chain_prefix = "Ch";
    total_prefix = "Tot";
  } else {
    chain_prefix = "Chain";
    total_prefix = "Total";
  }
}

void ProgressManager::maybePadToLength(std::string& content) const {
  if (!isRStudio) return;

  // Pad each line to exactly lineWidth characters (before adding \n)
  if (content.length() < lineWidth) {
    content += std::string(lineWidth - content.length(), ' ');
  } else if (content.length() > lineWidth) {
    content = content.substr(0, lineWidth);
  }
}


// Example usage/ test with RcppParallel
#include <RcppParallel.h>
// Worker functor for RcppParallel
struct ChainWorker : public RcppParallel::Worker {
  int nIter;
  ProgressManager &pm;
  bool display_progress;

  ChainWorker(int nIter_, ProgressManager &pm_, bool display_progress_)
    : nIter(nIter_), pm(pm_), display_progress(display_progress_) {}

  void operator()(std::size_t begin, std::size_t end) {

    auto chainId = begin;

    for (int i = 0; i < nIter; i++) {
      // ---- Simulated work ----
      std::this_thread::sleep_for(std::chrono::milliseconds(20));

      // ---- Update state ----
      pm.update(chainId);
      if (pm.shouldExit()) break;
      // if (Progress::check_abort()) Rcpp::checkUserInterrupt();
    }
  }
};


// [[Rcpp::export]]
void runMCMC_parallel(int nChains = 4, int nIter = 100, int nWarmup = 100, int progress_type = 2, bool useUnicode = false) {

  int nTotal = nIter + nWarmup;
  ProgressManager pm(nChains, nTotal, nWarmup, 10, progress_type, useUnicode);
  ChainWorker worker(nTotal, pm, true);

  // Run each chain in parallel
  RcppParallel::parallelFor(0, nChains, worker);

  if (pm.shouldExit()) {
    Rcpp::Rcout << "\nComputation interrupted by user.\n";
  } else {
    Rcpp::Rcout << "\nAll chains finished!\n";
  }
}

