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

    ProgressManager(int nChains_, int nIter_, int nWarmup_, int printEvery_ = 10, int progress_type = 2, bool useUnicode_ = false);
    void update(int chainId);
    void finish();
    bool shouldExit() const;

private:

    void checkConsoleWidthChange();
    int getConsoleWidth();
    std::string formatProgressBar(int chainId, int current, int total, double fraction, bool isTotal = false);
    std::string formatTimeInfo(int elapsed, int eta);
    void setupTheme();
    size_t getVisualLength(const std::string& str);

    bool isWarmupPhase() const {
        for (auto c : progress)
            if (c < nWarmup)
                return true;
        return false;
    }
    bool isWarmupPhase(const int chain_id) const {
      return progress[chain_id] < nWarmup;
    }

    void print();

    void update_prefixes(int width);

    void maybePadToLength(std::string& content) const;

    // Configuration parameters
    int nChains;                    ///< Number of parallel chains
    int nIter;                      ///< TOTAL Iterations per chain
    int nWarmup;                    ///< Warmup iterations per chain
    int printEvery;                 ///< Print frequency
    int no_spaces_for_total;        ///< Spacing for total line alignment
    int lastPrintedLines = 0;       ///< Lines printed in last update
    int lastPrintedChars = 0;       ///< Characters printed in last update (RStudio)
    int consoleWidth = 80;          ///< Current console width
    int lineWidth = 80;             ///< Target line width for content
    int prevConsoleWidth = -1;      ///< Previous console width for change detection

    // Environment and state flags
    bool isRStudio = false;         ///< Whether running in RStudio console
    bool needsToExit = false;       ///< User interrupt flag
    bool widthChanged = false;      ///< Console width changed flag

    // Visual configuration
    int barWidth = 40;              ///< Progress bar width in characters
    int progress_type = 2;          ///< Progress bar style type (0 = "none", 1 = "total", 2 = "per-chain")
    bool useUnicode = true;         ///< Use Unicode vs ASCII theme

    // Theme tokens
    std::string lhsToken;           ///< Left bracket/delimiter
    std::string rhsToken;           ///< Right bracket/delimiter
    std::string filledToken;        ///< Filled progress character
    std::string emptyToken;         ///< Empty progress character
    std::string partialTokenMore;   ///< Partial progress (>50%)
    std::string partialTokenLess;   ///< Partial progress (<50%)
    std::string chain_prefix;       ///< Chain label prefix
    std::string total_prefix;       ///< Total label prefix
    std::string total_padding;      ///< Padding for total line alignment

    // progress tracking
    std::vector<int> progress; ///< Per-chain progress counters

    // Timing
    Clock::time_point start;                                ///< Start time
    std::chrono::time_point<Clock> lastPrint;  ///< Last print time

    // Thread synchronization
    std::mutex printMutex;          ///< Mutex for thread-safe printing
};

#endif // PROGRESS_MANAGER_H