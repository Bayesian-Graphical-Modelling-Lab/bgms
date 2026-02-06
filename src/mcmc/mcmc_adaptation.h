// mcmc_adaptation.h
#pragma once

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "mcmc/mcmc_utils.h"
#include "mcmc/mcmc_rwm.h"
#include "math/explog_switch.h"

class DualAveraging {
public:
  double log_step_size;
  double log_step_size_avg;
  double hbar;
  double mu;
  double gamma;
  double t0;
  double kappa;
  int t;

  DualAveraging(double initial_step_size)
    : log_step_size(MY_LOG(initial_step_size)),
      log_step_size_avg(MY_LOG(initial_step_size)),
      hbar(0.0),
      mu(MY_LOG(10.0 * initial_step_size)),
      gamma(0.05),
      t0(10.0),
      kappa(0.75),
      t(1) {}

  void update(double accept_prob, double target_accept) {
    double eta = 1.0 / (t + t0);
    double error = target_accept - accept_prob;
    hbar = (1 - eta) * hbar + eta * error;
    log_step_size = mu - std::sqrt(t) / gamma * hbar;

    double weight = std::pow(t, -kappa);
    log_step_size_avg = weight * log_step_size + (1.0 - weight) * log_step_size_avg;
    t++;
  }

  void restart(double new_step_size) {
    log_step_size = MY_LOG(new_step_size);
    log_step_size_avg = MY_LOG(new_step_size);
    mu = MY_LOG(10.0 * new_step_size);
    hbar = 0.0;
    t = 1;
  }

  double current() const { return MY_EXP(log_step_size); }
  double averaged() const { return MY_EXP(log_step_size_avg); }
};



// === Diagonal Mass Matrix Estimator ===
class DiagMassMatrixAccumulator {
public:
  int count;
  arma::vec mean;
  arma::vec m2;

  DiagMassMatrixAccumulator(int dim)
    : count(0), mean(arma::zeros(dim)), m2(arma::zeros(dim)) {}

  void update(const arma::vec& sample) {
    count++;
    arma::vec delta = sample - mean;
    mean += delta / count;
    arma::vec delta2 = sample - mean;
    m2 += delta % delta2;
  }

  arma::vec variance() const {
    static constexpr double prior_weight = 5.0;
    static constexpr double prior_variance = 1e-3;
    double n = static_cast<double>(count);

    arma::vec empirical = m2 / std::max(1.0, n - 1.0);
    arma::vec prior = arma::ones(empirical.n_elem) * prior_variance;
    arma::vec var = (n / (n + prior_weight)) * empirical
    + (prior_weight / (n + prior_weight)) * prior;
    return var;
  }

  void reset() {
    count = 0;
    mean.zeros();
    m2.zeros();
  }
};


// === Dynamic Warmup Schedule with Adaptive Windows ===
//
// For edge_selection = FALSE:
//   Stage 1 (init), Stage 2 (doubling windows), Stage 3a (terminal)
//   total_warmup = user-specified warmup
//
// For edge_selection = TRUE:
//   User warmup is split: 85% for Stage 1-3a, 10% for Stage 3b, 5% for Stage 3c
//   - Stage 3b: proposal SD tuning for edge selection MH moves
//   - Stage 3c: step size re-adaptation with edge selection active
//   If Stage 3b would get < 20 iterations, it's skipped (uses default proposal SD)
//
// Warning types:
//   0 = none
//   1 = warmup extremely short (< 50)
//   2 = core stages using proportional fallback
//   3 = limited proposal SD tuning (edge_selection && warmup < 300)
//   4 = Stage 3b skipped (would have < 20 iterations)
//
struct WarmupSchedule {
  /* ---------- members ---------- */
  int stage1_end;                 // Stage-1   [0 … stage1_end-1]
  std::vector<int> window_ends;   // Stage-2   windows (last index of each)
  int stage3a_start;              // first iter in Stage-3 a
  int stage3b_start;              // first iter in Stage-3 b (== stage3c_start if skipped)
  int stage3c_start;              // first iter in Stage-3 c (== total_warmup if skipped)
  int total_warmup;               // warm-up iterations = user-specified value
  bool learn_proposal_sd;         // do we run the proposal-SD tuner?
  bool enable_selection;          // allow edge-indicator moves
  int warning_type;               // see above
  bool stage3b_skipped;           // true if 3b was skipped due to insufficient budget
  /* ------------------------------ */

  WarmupSchedule(int warmup,
                 bool enable_sel,
                 bool learn_sd)
    : stage1_end(0)
    , window_ends()
    , stage3a_start(0)
    , stage3b_start(0)
    , stage3c_start(0)
    , total_warmup(warmup)        // User gets exactly what they specify
    , learn_proposal_sd(learn_sd)
    , enable_selection(enable_sel)
    , warning_type(0)
    , stage3b_skipped(false)
  {
    // ===== Step 1: Determine budget allocation =====
    int warmup_core;    // Budget for Stage 1-3a (mass matrix + step size)
    int stage3b_budget = 0;

    if (enable_sel && learn_sd) {
      // For edge selection models: split warmup as 85%/10%/5%
      warmup_core    = static_cast<int>(0.85 * warmup);
      stage3b_budget = static_cast<int>(0.10 * warmup);
      // Stage 3c gets remainder: warmup - warmup_core - stage3b_budget (~5%)

      // Check if Stage 3b has enough iterations to be meaningful
      if (stage3b_budget < 20) {
        // Skip Stage 3b entirely - will use default proposal SD
        stage3b_skipped = true;
        warning_type = 4;
        warmup_core = warmup;  // Give all warmup to core stages
        stage3b_budget = 0;
      } else if (warmup < 300) {
        // Marginal but runs - warn about limited tuning
        warning_type = 3;
      }
    } else {
      // No edge selection: all warmup goes to core stages 1-3a
      warmup_core = warmup;
    }

    // ===== Step 2: Set up core stages (1, 2, 3a) =====
    constexpr int default_init_buffer = 75;   // Stage-1: initial fast adaptation
    constexpr int default_term_buffer = 50;   // Stage-3a: final fast adaptation
    constexpr int default_base_window = 25;   // Stage-2: initial window size

    int init_buffer, term_buffer, base_window;

    if (warmup_core < 20) {
      // Too short for any meaningful adaptation
      if (warning_type == 0) warning_type = 1;  // Don't overwrite more specific warnings
      init_buffer = warmup_core;
      term_buffer = 0;
      base_window = 0;
    } else if (default_init_buffer + default_base_window + default_term_buffer > warmup_core) {
      // Not enough room for fixed buffers; fall back to proportional (15%/75%/10%)
      if (warning_type == 0) warning_type = 2;
      init_buffer = static_cast<int>(0.15 * warmup_core);
      term_buffer = static_cast<int>(0.10 * warmup_core);
      base_window = warmup_core - init_buffer - term_buffer;
    } else {
      // Standard case: use fixed buffers
      init_buffer = default_init_buffer;
      term_buffer = default_term_buffer;
      base_window = default_base_window;
    }

    // Additional warning for extremely short warmup with edge selection
    if (enable_sel && warmup < 50 && warning_type != 1) {
      warning_type = 1;  // Override to most severe
    }

    /* ---------- Stage-1 ---------- */
    stage1_end = init_buffer;

    /* ---------- Stage-3a start ---------- */
    stage3a_start = warmup_core - term_buffer;

    /* ---------- Stage-2: build doubling windows ---------- */
    if (base_window > 0 && stage3a_start > stage1_end) {
      int cur   = stage1_end;
      int wsize = base_window;
      while (cur < stage3a_start) {
        int win = std::min(wsize, stage3a_start - cur);
        window_ends.push_back(cur + win);
        cur   += win;
        wsize  = std::min(wsize * 2, stage3a_start - cur);
      }
    }

    /* ---------- Stage-3b and 3c boundaries ---------- */
    stage3b_start = warmup_core;
    stage3c_start = warmup_core + stage3b_budget;
    // total_warmup already set to user's warmup value
  }

  /* ---------- helpers ---------- */
  bool in_stage1 (int i) const { return i <  stage1_end; }
  bool in_stage2 (int i) const { return i >= stage1_end && i < stage3a_start; }
  bool in_stage3a(int i) const { return i >= stage3a_start && i < stage3b_start; }
  bool in_stage3b(int i) const { return !stage3b_skipped && i >= stage3b_start && i < stage3c_start; }
  bool in_stage3c(int i) const { return enable_selection && !stage3b_skipped && i >= stage3c_start && i < total_warmup; }
  bool sampling (int i) const { return i >= total_warmup; }

  bool has_warning() const { return warning_type > 0; }
  bool warmup_extremely_short() const { return warning_type == 1; }
  bool using_proportional_fallback() const { return warning_type == 2; }
  bool limited_proposal_tuning() const { return warning_type == 3; }
  bool proposal_tuning_skipped() const { return warning_type == 4 || stage3b_skipped; }

  /* indicator moves in Stage 3c and sampling (if selection is enabled) */
  bool selection_enabled(int i) const {
    return enable_selection && (in_stage3c(i) || sampling(i));
  }

  /* adapt proposal_sd only inside Stage-3b (if that phase exists and wasn't skipped) */
  bool adapt_proposal_sd(int i) const {
    return learn_proposal_sd && !stage3b_skipped && in_stage3b(i);
  }

  /* current Stage-2 window index (-1 outside Stage-2) */
  int current_window(int i) const {
    for (size_t k = 0; k < window_ends.size(); ++k)
      if (i < window_ends[k]) return static_cast<int>(k);
    return -1;
  }
};



// === HMC/NUTS Adaptation Controller ===
class HMCAdaptationController {
public:
  HMCAdaptationController(int dim,
                          double initial_step_size,
                          double target_accept,
                          WarmupSchedule& schedule_ref,
                          bool learn_mass_matrix = true)
    : schedule(schedule_ref),
      learn_mass_matrix_(learn_mass_matrix),
      mass_accumulator(dim),
      step_adapter(initial_step_size),
      inv_mass_(arma::ones<arma::vec>(dim)),
      step_size_(initial_step_size),
      target_accept_(target_accept),
      finalized_mass_(false),
      mass_matrix_updated_(false) {}

  void update(const arma::vec& theta,
              double accept_prob,
              int iteration) {
    /* ---------------------------------------------------------
     * 1. STEP-SIZE ADAPTATION
     *    – runs in Stage-1, Stage-2, Stage-3a, and Stage-3c
     *    – Stage-3c continues adaptation because selection changes
     *      the model structure (including/excluding parameters)
     *    – pauses in Stage-3b (proposal SD learning phase)
     * --------------------------------------------------------- */
    if (schedule.in_stage1(iteration)  || schedule.in_stage2(iteration)  ||
    schedule.in_stage3a(iteration) || schedule.in_stage3c(iteration) )
    {
      step_adapter.update(accept_prob, target_accept_);
      step_size_ = step_adapter.current();
    }

    /* ---------------------------------------------------------
     * 2. MASS-MATRIX ADAPTATION
     *    – only while we are inside Stage-2
     * --------------------------------------------------------- */
    if (schedule.in_stage2(iteration) && learn_mass_matrix_) {
      mass_accumulator.update(theta);
      int w = schedule.current_window(iteration);
      if (iteration + 1 == schedule.window_ends[w]) {
        // STAN convention: inv_mass = variance (not 1/variance!)
        // Higher variance → higher inverse mass → parameter moves more freely
        inv_mass_ = mass_accumulator.variance();
        mass_accumulator.reset();
        // Signal that mass matrix was updated - caller should run heuristic
        // and call reinit_stepsize() with the new step size
        mass_matrix_updated_ = true;
      }
    }

    /* ---------------------------------------------------------
     * 3. FREEZE ε AS SOON AS WE ENTER STAGE-3b or SAMPLING
     * --------------------------------------------------------- */
    if (iteration == schedule.stage3b_start || schedule.sampling(iteration)) {
      step_size_ = step_adapter.averaged();
    }
  }

  double current_step_size() const { return step_size_; }
  double final_step_size() const { return step_adapter.averaged(); }
  const arma::vec& inv_mass_diag() const { return inv_mass_; }
  bool has_fixed_mass_matrix() const { return finalized_mass_; }

  /**
   * Check if the mass matrix was just updated and needs step size re-initialization.
   * After calling this, call reinit_stepsize() with the result of the heuristic.
   */
  bool mass_matrix_just_updated() const { return mass_matrix_updated_; }

  /**
   * Reinitialize step size adaptation after mass matrix update.
   * This should be called after running heuristic_initial_step_size() with
   * the new mass matrix to find an appropriate starting step size.
   *
   * - Set the new step size
   * - Set mu = log(10 * new_step_size) as the adaptation target
   * - Restart the dual averaging counters
   */
  void reinit_stepsize(double new_step_size) {
    step_size_ = new_step_size;
    step_adapter.restart(new_step_size);
    // Set mu to log(10 * epsilon) as per STAN's approach
    step_adapter.mu = MY_LOG(10.0 * new_step_size);
    mass_matrix_updated_ = false;
  }

private:
  WarmupSchedule& schedule;
  bool learn_mass_matrix_;
  DiagMassMatrixAccumulator mass_accumulator;
  DualAveraging step_adapter;
  arma::vec inv_mass_;
  double step_size_;
  double target_accept_;
  bool finalized_mass_;
  bool mass_matrix_updated_;
};



// RWMAdaptationController
// Handles Robbins-Monro adaptation for random walk Metropolis proposals
// Automatically adjusts step sizes during warmup using acceptance probabilities
class RWMAdaptationController {
public:
  arma::mat& proposal_sd;
  const int total_warmup;
  const double target_accept;

  RWMAdaptationController(arma::mat& proposal_sd_matrix,
                          const WarmupSchedule& warmup,
                          double target_accept_rate = 0.44)
    : proposal_sd(proposal_sd_matrix),
      total_warmup(warmup.total_warmup),
      target_accept(target_accept_rate) {}

  void update(const arma::umat& index_mask,
              const arma::mat& accept_prob_matrix,
              int iteration) {

    if (iteration >= total_warmup || iteration < 1)
      return;

    const double rm_decay_rate = 0.75;
    const double rm_weight = std::pow(iteration, -rm_decay_rate);

    for (arma::uword i = 0; i < proposal_sd.n_rows; ++i) {
      for (arma::uword j = 0; j < proposal_sd.n_cols; ++j) {
        if (index_mask(i, j) == 1) {
          const double accept_log = MY_LOG(accept_prob_matrix(i, j));
          proposal_sd(i, j) = update_proposal_sd_with_robbins_monro(
            proposal_sd(i, j),
            accept_log,
            rm_weight,
            target_accept
          );
        }
      }
    }
  }
};