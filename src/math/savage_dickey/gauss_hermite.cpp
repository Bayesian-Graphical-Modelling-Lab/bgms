#include "gauss_hermite.h"

#include <cmath>
#include <mutex>
#include <unordered_map>

namespace savage_dickey {

namespace {

// Precomputed Gauss-Hermite nodes/weights for several N (physicists',
// weight e^{-y²}). Built lazily via Golub-Welsch on first request and cached
// in a thread-safe static map, so per-call cost is amortised to a hashed
// lookup + N kernel evaluations.
struct GHRule {
    arma::vec nodes;
    arma::vec weights;
};

GHRule build_rule(int N) {
    arma::mat J(N, N, arma::fill::zeros);
    for (int k = 1; k < N; ++k) {
        const double off = std::sqrt(static_cast<double>(k) / 2.0);
        J(k - 1, k) = off;
        J(k, k - 1) = off;
    }
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, J);
    GHRule r;
    r.nodes   = eigval;
    r.weights = std::sqrt(arma::datum::pi) * arma::square(eigvec.row(0).t());
    return r;
}

const GHRule& get_rule(int N) {
    static std::unordered_map<int, GHRule> cache;
    static std::mutex mtx;
    std::lock_guard<std::mutex> lock(mtx);
    auto it = cache.find(N);
    if (it == cache.end()) {
        auto [ins, ok] = cache.emplace(N, build_rule(N));
        it = ins;
    }
    return it->second;
}

inline double log_kernel(double x, double A, double B,
                         double s, double alpha) {
    const double a1 = alpha - 1.0;
    double v = -A * x * x + B * x;
    if (a1 != 0.0) v += a1 * std::log(s + x * x);
    return v;
}

}  // namespace

LSDQuadResult density_at_l_ji_gh(double x_eval, double A, double B,
                                  double s_jj, double alpha,
                                  int num_nodes) {
    LSDQuadResult out;
    out.log_density = arma::datum::nan;
    out.log_Z       = arma::datum::nan;
    out.status      = 0;
    if (!(A > 0.0)) { out.status = 1; return out; }

    const GHRule& gh = get_rule(num_nodes);
    const double inv_sqrtA = 1.0 / std::sqrt(A);
    const double offset    = B / (2.0 * A);                  // Laplace center
    const double a1        = alpha - 1.0;

    // Compute log Σ_k w_k (s + (y_k/√A + B/(2A))²)^{a1}
    // via log-sum-exp on log(w_k) + a1·log(s + (...)²).
    arma::vec logterms(num_nodes);
    for (int k = 0; k < num_nodes; ++k) {
        const double x_k     = gh.nodes(k) * inv_sqrtA + offset;
        const double s_plus  = s_jj + x_k * x_k;
        logterms(k) = std::log(gh.weights(k))
                    + (a1 == 0.0 ? 0.0 : a1 * std::log(s_plus));
    }
    const double max_lt = logterms.max();
    const double lse = max_lt + std::log(arma::accu(arma::exp(logterms - max_lt)));

    // log Z = -½ log A + B²/(4A) + lse
    out.log_Z       = -0.5 * std::log(A) + (B * B) / (4.0 * A) + lse;
    out.log_density = log_kernel(x_eval, A, B, s_jj, alpha) - out.log_Z;
    out.status      = 0;
    return out;
}

}  // namespace savage_dickey
