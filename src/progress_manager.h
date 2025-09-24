#ifndef PROGRESS_MANAGER_H
#define PROGRESS_MANAGER_H

#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <numeric>

using Clock = std::chrono::steady_clock;

// Interrupt checking functions
// https://github.com/kforner/rcpp_progress/blob/d851ac62fd0314239e852392de7face5fa4bf48e/inst/include/interrupts.hpp#L24-L31
static void chkIntFn(void *dummy) {
	R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
inline bool checkInterrupt() {
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

/**
 * @brief Multi-chain progress bar manager for MCMC computations
 *
 * This class provides a thread-safe progress bar that works in both RStudio
 * console and terminal environments. It supports Unicode theming with colored
 * progress indicators and proper cursor positioning.
 *
 * Key features:
 * - Multi-chain progress tracking with atomic operations
 * - RStudio vs terminal environment detection and adaptation
 * - Unicode and classic theming options
 * - ANSI color support with proper visual length calculations
 * - Thread-safe printing with mutex protection
 * - Console width adaptation and change detection
 * - User interrupt checking
 */
class ProgressManager {

public:

    ProgressManager(int nChains_, int nIter_, int nWarmup_, int printEvery_ = 10, int progress_type = 2, bool useUnicode_ = true);
    void update(size_t chainId);
    void finish();
    bool shouldExit() const;

private:

    void checkConsoleWidthChange();
    size_t getConsoleWidth() const;
    std::string formatProgressBar(size_t chainId, size_t current, size_t total, double fraction, bool isTotal = false) const;
    std::string formatTimeInfo(double elapsed, double eta) const;
    std::string formatDuration(double seconds) const;
    void setupTheme();

    bool isWarmupPhase() const {
        for (auto c : progress)
            if (c < nWarmup)
                return true;
        return false;
    }
    bool isWarmupPhase(const size_t chain_id) const {
      return progress[chain_id] < nWarmup;
    }

    void print();

    void update_prefixes(size_t  width);

    void maybePadToLength(std::string& content) const;

    // set by constructor
    size_t nChains;                    // Number of parallel chains
    size_t nIter;                      // Total Iterations per chain
    size_t nWarmup;                    // Warmup iterations per chain
    size_t printEvery;                 // Print frequency
    size_t progress_type = 2;          // Progress bar style type (0 = "none", 1 = "total", 2 = "per-chain")
    bool useUnicode = true;            // Use Unicode vs ASCII theme
    std::vector<size_t> progress;      // Per-chain progress counters

    // internal config parameters/ data
    size_t no_spaces_for_total;     // Spacing for total line alignment
    size_t lastPrintedLines = 0;    // Lines printed in last update
    size_t lastPrintedChars = 0;    // Characters printed in last update (RStudio)
    size_t consoleWidth = 80;       // Current console width
    size_t lineWidth = 80;          // Target line width for content
    int prevConsoleWidth = -1;      // Previous console width for change detection

    // Environment and state flags
    bool isRStudio = false;         ///< Whether running in RStudio console
    bool needsToExit = false;       ///< User interrupt flag
    bool widthChanged = false;      ///< Console width changed flag

    // Visual configuration
    size_t barWidth = 40;              // Progress bar width in characters

    // Theme tokens
    std::string lhsToken;           // Left bracket/delimiter
    std::string rhsToken;           // Right bracket/delimiter
    std::string filledToken;        // Filled progress character
    std::string emptyToken;         // Empty progress character
    std::string partialTokenMore;   // Partial progress (>50%)
    std::string partialTokenLess;   // Partial progress (<50%)
    std::string chain_prefix;       // Chain label prefix
    std::string total_prefix;       // Total label prefix
    std::string total_padding;      // Padding for total line alignment

    // Timing
    Clock::time_point start;                   // Start time
    std::chrono::time_point<Clock> lastPrint;  // Last print time

    // Thread synchronization
    std::mutex printMutex;          // Mutex for thread-safe printing
};

#endif // PROGRESS_MANAGER_H