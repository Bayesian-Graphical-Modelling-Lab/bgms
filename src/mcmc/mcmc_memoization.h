#pragma once

#include <RcppArmadillo.h>
#include <functional>


/**
 * Class: Memoizer
 *
 * Single-entry cache for joint log-posterior and gradient evaluations.
 *
 * In NUTS, the typical access pattern within a leapfrog step is:
 *   1. cached_grad(theta) — compute gradient (and cache logp as side-effect)
 *   2. cached_log_post(theta) — retrieve the already-cached logp
 *
 * A single-entry cache is optimal here because each leapfrog step produces
 * a new unique theta: hash-map lookups would almost never hit, and hashing
 * an arma::vec element-by-element is expensive.
 *
 * The joint evaluation function computes both logp and gradient together
 * (since models often share most of the computation between the two).
 */
class Memoizer {
public:
  using JointFn = std::function<std::pair<double, arma::vec>(const arma::vec&)>;

  JointFn joint_fn;

  // Single-entry cache
  arma::vec cached_theta;
  double    cached_logp_val;
  arma::vec cached_grad_val;
  bool      has_cache = false;

  /**
   * Construct from separate log_post and grad functions.
   * Calls them independently (backward-compatible).
   */
  Memoizer(
    const std::function<double(const arma::vec&)>& lp,
    const std::function<arma::vec(const arma::vec&)>& gr
  ) : joint_fn([lp, gr](const arma::vec& theta) -> std::pair<double, arma::vec> {
        arma::vec g = gr(theta);
        double v = lp(theta);
        return {v, std::move(g)};
      }) {}

  /**
   * Construct from a joint function that computes both at once.
   */
  explicit Memoizer(JointFn jf) : joint_fn(std::move(jf)) {}

  double cached_log_post(const arma::vec& theta) {
    ensure_cached(theta);
    return cached_logp_val;
  }

  const arma::vec& cached_grad(const arma::vec& theta) {
    ensure_cached(theta);
    return cached_grad_val;
  }

private:
  void ensure_cached(const arma::vec& theta) {
    if (has_cache &&
        theta.n_elem == cached_theta.n_elem &&
        std::memcmp(theta.memptr(), cached_theta.memptr(),
                    theta.n_elem * sizeof(double)) == 0) {
      return;
    }
    auto [lp, gr] = joint_fn(theta);
    cached_theta = theta;
    cached_logp_val = lp;
    cached_grad_val = std::move(gr);
    has_cache = true;
  }
};