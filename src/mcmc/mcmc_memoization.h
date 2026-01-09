#pragma once

#include <RcppArmadillo.h>
#include <unordered_map>     // for std::unordered_map
#include <functional>        // for std::function



/**
 * Struct: VecHash
 *
 * A hash function is a way of converting complex input data (like a vector of numbers)
 * into a single fixed-size integer value. This hash acts like a digital fingerprint —
 * it helps uniquely identify the input and allows fast lookup in data structures
 * like hash tables.
 *
 * In C++, hash tables are implemented using containers like std::unordered_map,
 * which require a hash function and an equality comparator for the key type.
 * Since arma::vec (a vector of doubles) isn't natively supported as a key,
 * this struct defines how to compute a hash for it.
 *
 * Why use hashing here?
 *   - We cache (memoize) values like log-posterior or gradient at specific points.
 *   - To retrieve them quickly later, we need a fast way to index by vector.
 *   - Hashing gives constant-time access (on average) by turning the vector into a hash code.
 *
 * This implementation:
 *   - Initializes a seed based on the vector size.
 *   - Iteratively incorporates each element's hash into the seed using bitwise operations
 *     and a constant based on the golden ratio to reduce collisions.
 *
 * The result is a unique and repeatable hash value that lets us efficiently
 * store and retrieve cached evaluations.
 *
 * Provides a custom hash function for arma::vec so that it can be used as a key
 * in standard hash-based containers like std::unordered_map.
 *
 * Since arma::vec is not hashable by default in C++, this functor allows you
 * to compute a combined hash from the individual elements of the vector.
 *
 * Hashing Strategy:
 *   - Starts with an initial seed equal to the vector's length.
 *   - Iteratively combines the hash of each element into the seed using
 *     a standard hashing recipe (XOR, bit-shift mix, and golden ratio).
 *
 * This ensures that vectors with similar but not identical elements
 * produce distinct hash codes, reducing the risk of collisions.
 */
struct VecHash {
  std::size_t operator()(const arma::vec& v) const {
    std::size_t seed = v.n_elem;
    for (size_t i = 0; i < v.n_elem; i++) {
      seed ^= std::hash<double>{}(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};



/**
 * Struct: VecEqual
 *
 * Provides a custom equality comparator for arma::vec, to be used in conjunction
 * with VecHash in unordered_map and similar containers.
 *
 * By default, arma::vec does not implement equality suitable for use as a key
 * in hashed containers. This comparator solves that by relying on Armadillo's
 * approx_equal function.
 *
 * Comparison Strategy:
 *   - Uses "absdiff" mode to check if the absolute difference between vectors
 *     is below a specified tolerance.
 *   - The tolerance is hardcoded as 1e-12, assuming numerical precision errors
 *     should not cause false inequality.
 *
 * This is crucial for numerical stability, especially when the same vector
 * might be computed multiple times with slight rounding differences.
 */
struct VecEqual {
  bool operator()(const arma::vec& a, const arma::vec& b) const {
    return arma::approx_equal(a, b, "absdiff", 1e-12);
  }
};



/**
 * Class: Memoizer
 *
 * A utility class that caches (memoizes) evaluations of the log-posterior and its gradient
 * at specific parameter values (theta) to avoid redundant and costly recomputation.
 *
 * Usage:
 *   - Constructed with a log-posterior function and its gradient.
 *   - Uses hash-based caching (via unordered_map) keyed by arma::vec inputs.
 *     That means:
 *       - Each parameter vector (theta) is used as a key.
 *       - Results are stored in a lookup table so that repeated evaluations can be skipped.
 *       - A custom VecHash function generates a hash code for arma::vec.
 *       - A custom VecEqual comparator defines when two vectors are considered equal (with a small tolerance).
 *     This enables fast retrieval (average constant time) of previously computed values for log_post and grad.
 *
 * Methods:
 *  - cached_log_post(const arma::vec& theta):
 *        Returns the cached log-posterior if available, otherwise computes, stores, and returns it.
 *
 *  - cached_grad(const arma::vec& theta):
 *        Returns the cached gradient if available, otherwise computes, stores, and returns it.
 *
 * Internal Details:
 *  - Relies on VecHash and VecEqual to use arma::vec as map keys.
 *  - Assumes high-precision comparison for theta keys (absdiff tolerance 1e-12).
 */
class Memoizer {
public:
  std::function<double(const arma::vec&)> log_post;
  std::function<arma::vec(const arma::vec&)> grad;

  std::unordered_map<arma::vec, double, VecHash, VecEqual> logp_cache;
  std::unordered_map<arma::vec, arma::vec, VecHash, VecEqual> grad_cache;

  Memoizer(
    const std::function<double(const arma::vec&)>& lp,
    const std::function<arma::vec(const arma::vec&)>& gr
  ) : log_post(lp), grad(gr) {}

  double cached_log_post(const arma::vec& theta) {
    auto it = logp_cache.find(theta);
    if (it != logp_cache.end()) return it->second;
    double val = log_post(theta);
    logp_cache[theta] = val;
    return val;
  }

  arma::vec cached_grad(const arma::vec& theta) {
    auto it = grad_cache.find(theta);
    if (it != grad_cache.end()) return it->second;
    arma::vec val = grad(theta);
    grad_cache[theta] = val;
    return val;
  }
};