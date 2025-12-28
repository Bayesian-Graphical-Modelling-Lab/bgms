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

  arma::vec inverse_mass() const {
    return 1.0 / variance();
  }

  void reset() {
    count = 0;
    mean.zeros();
    m2.zeros();
  }
};


// === Stan-style Dynamic Warmup Schedule with Adaptive Windows ===
struct WarmupSchedule {
  /* ---------- members (keep this order!) ---------- */
  int stage1_end;                 // Stage-1   [0 … stage1_end-1]
  std::vector<int> window_ends;   // Stage-2   windows (last index of each)
  int stage3a_start;              // first iter in Stage-3 a
  int stage3b_start;              // first iter in Stage-3 b (may == total_warmup)
  int stage3c_start;              // = stage3b_end if selection is ON, else stage3b_end == total_warmup
  int total_warmup;               // warm-up iterations in total
  bool learn_proposal_sd;         // do we run the proposal-SD tuner?
  bool enable_selection;          // allow edge-indicator moves (after warm-up)
  /* ------------------------------------------------ */

  WarmupSchedule(int warmup_core,
                 bool enable_sel,
                 bool learn_sd)
    : stage1_end(0)
    , window_ends()
    , stage3a_start(0)
    , stage3b_start(0)
    , stage3c_start(0)
    , total_warmup(0)             // will be set below
    , learn_proposal_sd(learn_sd)
    , enable_selection(enable_sel)
  {
    /* ---------- Stage-1 (7.5 % of core warm-up) ---------- */
    stage1_end      = int(0.075 * warmup_core);

    /* ---------- Stage-3 a (last 10 % of core warm-up) ---- */
    stage3a_start   = warmup_core - int(0.10 * warmup_core);

    /* ---------- Stage-2 : build doubling windows ---------- */
    int cur   = stage1_end;
    int wsize = 25;
    while (cur < stage3a_start) {
      int win = std::min(wsize, stage3a_start - cur);
      window_ends.push_back(cur + win);          // 0-based index of *last* iter
      cur   += win;
      wsize  = std::min(wsize * 2, stage3a_start - cur);
    }

    /* ---------- Stage-3 b : proposal-SD learning ---------- */
    int stage3a_len = warmup_core - stage3a_start;          // 10 % of core
    int stage3b_len = (learn_proposal_sd && enable_selection)
      ? std::max(100, stage3a_len)
        : 0;
    int stage3c_len = (learn_proposal_sd && enable_selection)
      ? std::max(100, stage3a_len)
        : 0;

    stage3b_start = warmup_core;            // begins directly after 3 a
    stage3c_start = warmup_core + stage3b_len;
    total_warmup  = stage3c_start + stage3c_len;
  }

  /* ---------- helpers ---------- */
  bool in_stage1 (int i) const { return i <  stage1_end;               }
  bool in_stage2 (int i) const { return i >= stage1_end   && i < stage3a_start; }
  bool in_stage3a(int i) const { return i >= stage3a_start&& i < stage3b_start; }
  bool in_stage3b(int i) const { return i >= stage3b_start&& i < stage3c_start; }
  bool in_stage3c(int i) const { return enable_selection && i >= stage3c_start && i < total_warmup; }
  bool sampling (int i) const { return i >= total_warmup; }

  /* indicator moves only once warm-up is finished */
  bool selection_enabled(int i) const {
    return enable_selection && ( in_stage3c(i) || sampling(i) );
  }

  /* adapt proposal_sd only inside Stage-3 b (if that phase exists) */
  bool adapt_proposal_sd(int i) const {
    return learn_proposal_sd && in_stage3b(i);
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
      finalized_mass_(false) {}

  void update(const arma::vec& theta,
              double accept_prob,
              int iteration) {
    /* ---------------------------------------------------------
     * 1. STEP-SIZE ADAPTATION
     *    – should run in Stage-1, Stage-2, Stage-3a
     *    – must stop in Stage-3b and afterwards
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
        inv_mass_ = 1.0 / mass_accumulator.variance();
        mass_accumulator.reset();
        step_adapter.restart(step_size_);
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

private:
  WarmupSchedule& schedule;
  bool learn_mass_matrix_;
  DiagMassMatrixAccumulator mass_accumulator;
  DualAveraging step_adapter;
  arma::vec inv_mass_;
  double step_size_;
  double target_accept_;
  bool finalized_mass_;
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