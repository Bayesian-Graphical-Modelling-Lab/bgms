#include "sd_density_l_space_quad.h"
#include "sd_density_cubic.h"

#include <cmath>
#include <mutex>
#include <unordered_map>

namespace ggm_sd {

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


namespace {

// Genz-Keister nested Kronrod-Patterson rules for the physicists' Hermite
// weight exp(-y^2). Three strictly-positive-weight levels: 3, 9, 35 nodes.
// Nodes nest: GK_3 nodes are a subset of GK_9, which is a subset of GK_35.
// Tables extracted from R's SparseGrid (type "KPN", k = 2, 5, 18) and
// converted from probabilist weight via y = x_prob / sqrt(2),
// w_phys = sqrt(pi) * w_prob. Polynomial exactness (verified against even
// moments int y^(2m) exp(-y^2) dy = (2m-1)!! / 2^m sqrt(pi)):
//   GK_3:  degree 5
//   GK_9:  degree 15
//   GK_35: degree at least 50.

static constexpr int   kGK3Nodes   = 3;
static constexpr double kGK3Nodes_y[3] = {
    -1.22474487139158894e+00,
    0.00000000000000000e+00,
    1.22474487139158894e+00
};
static constexpr double kGK3Nodes_w[3] = {
    2.95408975150919351e-01,
    1.18163590060367740e+00,
    2.95408975150919351e-01
};

static constexpr int   kGK9Nodes   = 9;
static constexpr double kGK9Nodes_y[9] = {
    -2.95921077906383800e+00,
    -2.02323019110051572e+00,
    -1.22474487139158894e+00,
    -5.24033547486957629e-01,
    0.00000000000000000e+00,
    5.24033547486957629e-01,
    1.22474487139158894e+00,
    2.02323019110051572e+00,
    2.95921077906383800e+00
};
static constexpr double kGK9Nodes_w[9] = {
    1.67088263068823480e-04,
    1.41731178739790977e-02,
    1.68118928947677706e-01,
    4.78694285491141238e-01,
    4.50147009753781968e-01,
    4.78694285491141238e-01,
    1.68118928947677706e-01,
    1.41731178739790977e-02,
    1.67088263068823480e-04
};

static constexpr int   kGK35Nodes  = 35;
static constexpr double kGK35Nodes_y[35] = {
    -6.37593927098223561e+00,
    -5.64325785788574393e+00,
    -5.03608994447309401e+00,
    -4.49959939831038811e+00,
    -4.02922014050437129e+00,
    -3.66777421594633781e+00,
    -3.34916395371319446e+00,
    -2.95921077906383800e+00,
    -2.57055837658429683e+00,
    -2.26651326205678759e+00,
    -2.02323019110051572e+00,
    -1.83570797517518680e+00,
    -1.57941213484676712e+00,
    -1.22474487139158894e+00,
    -8.70040895352902852e-01,
    -5.24033547486957629e-01,
    -1.76064142082008934e-01,
    0.00000000000000000e+00,
    1.76064142082008934e-01,
    5.24033547486957629e-01,
    8.70040895352902852e-01,
    1.22474487139158894e+00,
    1.57941213484676712e+00,
    1.83570797517518680e+00,
    2.02323019110051572e+00,
    2.26651326205678759e+00,
    2.57055837658429683e+00,
    2.95921077906383800e+00,
    3.34916395371319446e+00,
    3.66777421594633781e+00,
    4.02922014050437129e+00,
    4.49959939831038811e+00,
    5.03608994447309401e+00,
    5.64325785788574393e+00,
    6.37593927098223561e+00
};
static constexpr double kGK35Nodes_w[35] = {
    1.86840148945106082e-18,
    9.65994662785632742e-15,
    5.48968369484994948e-12,
    8.15537218169169284e-10,
    3.79202223923195514e-08,
    4.37378180409270046e-07,
    4.84627997370204776e-06,
    6.33286208056179050e-05,
    4.87853993044437865e-04,
    1.45155804251559102e-03,
    4.09675277203440640e-03,
    5.59288289114691969e-03,
    2.77805089085351106e-02,
    8.02455181473909207e-02,
    1.63712215557358126e-01,
    2.62448714887842827e-01,
    3.39885955855852240e-01,
    9.12626753637379540e-04,
    3.39885955855852240e-01,
    2.62448714887842827e-01,
    1.63712215557358126e-01,
    8.02455181473909207e-02,
    2.77805089085351106e-02,
    5.59288289114691969e-03,
    4.09675277203440640e-03,
    1.45155804251559102e-03,
    4.87853993044437865e-04,
    6.33286208056179050e-05,
    4.84627997370204776e-06,
    4.37378180409270046e-07,
    3.79202223923195514e-08,
    8.15537218169169284e-10,
    5.48968369484994948e-12,
    9.65994662785632742e-15,
    1.86840148945106082e-18
};

// AGHQ log Z at a single Genz-Keister level. Centres the N nested nodes at
// (phi_star, kappa) and returns the log-sum-exp normaliser. The cascade
// driver calls this at increasing N. All weights are strictly positive so
// the standard log-sum-exp applies.
inline double aghq_log_Z_at_level(int num_nodes,
                                  const double* nodes_y,
                                  const double* nodes_w,
                                  double phi_star, double kappa,
                                  double A, double B, double s_jj,
                                  double alpha) {
    const double scale     = std::sqrt(2.0 / kappa);
    const double log_scale = 0.5 * std::log(2.0 / kappa);
    arma::vec logterms(num_nodes);
    for (int k = 0; k < num_nodes; ++k) {
        const double y_k   = nodes_y[k];
        const double phi_k = phi_star + y_k * scale;
        logterms(k) = std::log(nodes_w[k])
                    + log_kernel(phi_k, A, B, s_jj, alpha)
                    + y_k * y_k;
    }
    const double max_lt = logterms.max();
    const double lse = max_lt + std::log(arma::accu(arma::exp(logterms - max_lt)));
    return log_scale + lse;
}

}  // namespace

LSDAGHQResult density_at_l_ji_aghq(double x_eval, double A, double B,
                                    double s_jj, double alpha) {
    LSDAGHQResult out;
    out.log_density   = arma::datum::nan;
    out.log_Z         = arma::datum::nan;
    out.log_Z_err_est = arma::datum::nan;
    out.x_mode        = arma::datum::nan;
    out.curvature     = arma::datum::nan;
    out.n_nodes_used  = 0;
    out.status        = 0;
    if (!(A > 0.0))    { out.status = 1; return out; }
    if (!(s_jj > 0.0)) { out.status = 2; return out; }

    // Closed-form mode lookup via the critical-point cubic. Fall back to
    // the alpha = 1 reference Gaussian (B/(2A), 2A) when the solver returns
    // no usable mode (PD revert, triple-root degeneracy, or curvature below
    // kappa_floor = 1e-10 * 2A from roundoff that promotes ell_pp = 0 to
    // a slightly negative value).
    const CubicResult cubic = solve_sd_cubic(A, B, s_jj, alpha);
    const double kappa_floor = 1e-10 * 2.0 * A;
    double phi_star, kappa;
    bool   used_fallback = false;
    if (cubic.status == 0
        && cubic.n_modes >= 1
        && cubic.global_mode_index >= 0
        && cubic.roots[cubic.global_mode_index].curvature > kappa_floor) {
        phi_star = cubic.roots[cubic.global_mode_index].phi;
        kappa    = cubic.roots[cubic.global_mode_index].curvature;
    } else {
        phi_star      = B / (2.0 * A);
        kappa         = 2.0 * A;
        used_fallback = true;
    }
    out.x_mode    = phi_star;
    out.curvature = kappa;

    // Three-level Genz-Keister escalation cascade (nested, strictly
    // positive weights): GK_3 -> GK_9 -> GK_35. Polynomial exactness
    // rises 5 -> 15 -> 50+, so the cascade jumps through several orders
    // of magnitude of truncation error per step on analytic integrands.
    //
    // Tolerance choice: |I_high - I_low| is the *rule-difference* (the
    // cheaper rule's truncation error), not the higher rule's error
    // against truth. With GK_9 polynomial degree 15 and GK_35 degree
    // 50+, GK_35's actual error is many orders of magnitude smaller
    // than the GK_9-GK_35 rule difference. tol_strict = 1e-3 therefore
    // means "GK_9 was already converged to 1e-3" -- which by extension
    // means GK_35 is converged to roundoff. Cells where the cascade
    // fails to meet this threshold (rule_diff > 1e-3) flag genuine
    // non-Gaussian structure (bimodal, sharp log feature) where the
    // single-Laplace AGHQ at the global mode is suspect; the chain
    // PD-reverts on status = 4 in those cells.
    constexpr double kTolStrict = 1e-3;

    double I_prev = aghq_log_Z_at_level(kGK3Nodes, kGK3Nodes_y, kGK3Nodes_w,
                                         phi_star, kappa, A, B, s_jj, alpha);
    int    n_prev = kGK3Nodes;
    int    converged_level = -1;
    double err = arma::datum::nan;
    double I_final = I_prev;

    struct Level { int n; const double* y; const double* w; };
    constexpr int kNumEscalations = 2;
    const Level levels[kNumEscalations] = {
        {kGK9Nodes,  kGK9Nodes_y,  kGK9Nodes_w},
        {kGK35Nodes, kGK35Nodes_y, kGK35Nodes_w}
    };
    for (int lvl = 0; lvl < kNumEscalations; ++lvl) {
        const double I_cur = aghq_log_Z_at_level(
            levels[lvl].n, levels[lvl].y, levels[lvl].w,
            phi_star, kappa, A, B, s_jj, alpha);
        err     = std::abs(I_cur - I_prev);
        I_final = I_cur;
        n_prev  = levels[lvl].n;
        if (err < kTolStrict) {
            converged_level = levels[lvl].n;
            break;
        }
        I_prev = I_cur;
    }

    out.log_Z         = I_final;
    out.log_Z_err_est = err;
    out.n_nodes_used  = n_prev;
    out.log_density   = log_kernel(x_eval, A, B, s_jj, alpha) - out.log_Z;
    if (converged_level < 0) {
        // Cascade exhausted without converging within tol_strict; signal
        // status = 4 so the caller can PD-revert.
        out.status = 4;
    } else {
        out.status = used_fallback ? 3 : 0;
    }
    return out;
}

}  // namespace ggm_sd
