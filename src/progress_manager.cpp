#include "progress_manager.h"


// ProgressManager Implementation

ProgressManager::ProgressManager(int nChains_, int nIter_, int printEvery_, bool useUnicode_)
  : nChains(nChains_), nIter(nIter_), progress(nChains_), printEvery(printEvery_),
    lastLineLengths(nChains_ + 2, 0), useUnicode(useUnicode_) { // +2 for total and time lines

  for (int i = 0; i < nChains; i++) progress[i] = 0;
  start = Clock::now();
  lastPrint.store(start);

  // Check if we're in RStudio
  Environment base("package:base");
  Function getOption = base["getOption"];
  Function Sysgetenv = base["Sys.getenv"];
  SEXP s = Sysgetenv("RSTUDIO");
  isRStudio = as<std::string>(s) == "1";

  no_spaces_for_total = 3 + static_cast<int>(std::log10(nChains));

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

ProgressManager::~ProgressManager() {
  // Skip destructor print to avoid double printing
  // The final state is already shown by the last regular print
}

void ProgressManager::update(int chainId) {
  progress[chainId]++;
  totalDone++;

  // Only chain 0 handles printing to avoid race conditions
  if (chainId == 0 && progress[chainId] % printEvery == 0) {
    auto now = Clock::now();
    auto lastPrintTime = lastPrint.load();
    std::chrono::duration<double> sinceLast = now - lastPrintTime;

    // Throttle printing to avoid spamming
    if (sinceLast.count() >= 0.5) {
      print();
      lastPrint.store(now);
    }
  }

  // Check for user interrupts and console width changes less frequently to reduce overhead
  if (chainId == 0 && progress[chainId] % (printEvery * 5) == 0) {
    needsToExit = checkInterrupt();
    // Also check for console width changes occasionally
    checkConsoleWidthChange();
  }
}

bool ProgressManager::shouldExit() const {
  return needsToExit;
}

void ProgressManager::checkConsoleWidthChange() {
  if (!isRStudio) return;

  Environment base("package:base");
  Function getOption = base["getOption"];
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
  Environment base("package:base");
  Function getOption = base["getOption"];
  SEXP s = getOption("width", 0);
  int width = as<int>(s);
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
    builder << total_prefix << ":" << std::string(no_spaces_for_total, ' ') << progressBar << " " << current << "/" << total
            << " (" << std::fixed << std::setprecision(1) << fraction * 100 << "%)";
  } else {
      builder << chain_prefix << " " << chainId << ": " << progressBar << " " << current << "/" << total
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
    filledToken = "\x1b[34m━\x1b[0m";  // Blue filled
    emptyToken = "\x1b[37m━\x1b[0m";   // Gray empty
    partialTokenMore = "\x1b[34m╸\x1b[0m";  // Blue partial (> 0.5)
    partialTokenLess = "\x1b[37m╺\x1b[0m";  // Gray partial (< 0.5)
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
  int done = totalDone;
  double fracTotal = double(done) / totalWork;
  // should actually be the eta of the slowest chain!
  double eta = (fracTotal > 0) ? elapsed / fracTotal - elapsed : 0.0;

  std::ostringstream out;
  int totalChars = 0;
  int lineIndex = 0;

  if (isRStudio) {
    // if this is not the first print, delete previous content
    if (lastPrintedChars > 0) {
      out << "\x1b[" << std::to_string(lastPrintedChars) << "D";
    }

    std::string clearedLine;

    // Print progress for each chain
    for (int i = 0; i < nChains; i++) {
      double frac = double(progress[i]) / nIter;
      std::string chainProgress = formatProgressBar(i + 1, progress[i], nIter, frac);

      // Pad each line to exactly lineWidth characters (before adding \n)
      if (chainProgress.length() < lineWidth) {
        chainProgress += std::string(lineWidth - chainProgress.length(), ' ');
      } else if (chainProgress.length() > lineWidth) {
        chainProgress = chainProgress.substr(0, lineWidth);
      }

      // Add newline to every line (including the last one)
      clearedLine = chainProgress + "\n";
      out << clearedLine;

      lastLineLengths[lineIndex] = clearedLine.length(); // This will be lineWidth + 1
      totalChars += lastLineLengths[lineIndex];
      lineIndex++;
    }

    // Print total progress
    std::string totalProgress = formatProgressBar(0, done, totalWork, fracTotal, true);

    // Pad total progress line to exactly lineWidth characters (before adding \n)
    if (totalProgress.length() < lineWidth) {
      totalProgress += std::string(lineWidth - totalProgress.length(), ' ');
    } else if (totalProgress.length() > lineWidth) {
      totalProgress = totalProgress.substr(0, lineWidth);
    }

    clearedLine = totalProgress + "\n";
    out << clearedLine;

    lastLineLengths[lineIndex] = clearedLine.length(); // This will be lineWidth + 1
    totalChars += lastLineLengths[lineIndex];
    lineIndex++;

    // Print time info
    std::string timeInfo = formatTimeInfo(int(elapsed), int(eta));

    // Pad time info line to exactly lineWidth characters (before adding \n)
    if (timeInfo.length() < lineWidth) {
      timeInfo += std::string(lineWidth - timeInfo.length(), ' ');
    } else if (timeInfo.length() > lineWidth) {
      timeInfo = timeInfo.substr(0, lineWidth);
    }

    clearedLine = timeInfo + "\n";
    out << clearedLine;

    lastLineLengths[lineIndex] = clearedLine.length(); // This will be lineWidth + 1
    totalChars += lastLineLengths[lineIndex];

    // Track characters and lines for next cursor movement
    lastPrintedChars = totalChars;
    lastPrintedLines = lineIndex;

    std::string out_str = out.str();

    assert(out_str.length() == static_cast<size_t>(totalChars));
    Rcpp::Rcout << out_str;

  } else {
    // Terminal: Use carriage return and clear lines approach
    if (lastPrintedLines > 0) {
      // Move cursor up to start of our content and clear everything
      for (int i = 0; i < lastPrintedLines; i++) {
        out << "\x1b[1A\x1b[2K"; // Move up one line and clear entire line
      }
    }

    // Print progress for each chain
    for (int i = 0; i < nChains; i++) {
      double frac = double(progress[i]) / nIter;
      std::string chainProgress = formatProgressBar(i + 1, progress[i], nIter, frac);
      out << chainProgress << "\n";
    }

    // Print total progress
    std::string totalProgress = formatProgressBar(0, done, totalWork, fracTotal, true);
    out << totalProgress << "\n";

    // Print time info
    std::string timeInfo = formatTimeInfo(int(elapsed), int(eta));
    out << timeInfo << "\n";

    // Track total lines printed (chains + total + time)
    lastPrintedLines = nChains + 2;

    Rcpp::Rcout << out.str();
  }
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

std::string ProgressManager::addPadding(const std::string& content) {
  int paddingNeeded = consoleWidth - content.length();
  if (paddingNeeded > 0) {
    return content + std::string(paddingNeeded, ' ');
  } else {
    return content;
  }
}

std::string ProgressManager::clearLineLeftovers(const std::string& newContent, int oldLineLength) {
  if (newContent.length() < oldLineLength) {
    // Add spaces to clear old content
    return newContent + std::string(oldLineLength - newContent.length(), ' ');
  }
  return newContent;
}

// Worker functor for RcppParallel
struct ChainWorker : public Worker {
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
void runMCMC_parallel(int nChains = 4, int nIter = 100, bool display_progress = true, bool useUnicode = true) {

  ProgressManager pm(nChains, nIter, 10, useUnicode);
  ChainWorker worker(nIter, pm, display_progress);

  // Run each chain in parallel
  parallelFor(0, nChains, worker);

  if (pm.shouldExit()) {
    Rcpp::Rcout << "\nComputation interrupted by user.\n";
  } else {
    Rcpp::Rcout << "\nAll chains finished!\n";
  }
}

