#ifndef PROGRESS_MANAGER_H
#define PROGRESS_MANAGER_H

#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;
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

    ProgressManager(int nChains_, int nIter_, int printEvery_ = 10, bool useUnicode_ = false);
    ~ProgressManager();
    void update(int chainId);
    bool shouldExit() const;

private:

    void checkConsoleWidthChange();
    int getConsoleWidth();
    std::string formatProgressBar(int chainId, int current, int total, double fraction, bool isTotal = false);
    std::string formatTimeInfo(int elapsed, int eta);
    void setupTheme();
    size_t getVisualLength(const std::string& str);

    void print();

    void update_prefixes(int width);

    std::string addPadding(const std::string& content);

    std::string clearLineLeftovers(const std::string& newContent, int oldLineLength);

    // Configuration parameters
    int nChains;                    ///< Number of parallel chains
    int nIter;                      ///< Iterations per chain
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

    // Line tracking for cursor positioning
    std::vector<int> lastLineLengths; ///< Track length of each printed line

    // Thread-safe progress tracking
    std::vector<std::atomic<int>> progress; ///< Per-chain progress counters
    std::atomic<int> totalDone{0};          ///< Total completed iterations

    // Timing
    Clock::time_point start;                                ///< Start time
    std::atomic<std::chrono::time_point<Clock>> lastPrint;  ///< Last print time

    // Thread synchronization
    std::mutex printMutex;          ///< Mutex for thread-safe printing
};

#endif // PROGRESS_MANAGER_H