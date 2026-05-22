#pragma once

#include <memory>
#include <string>
#include <cmath>
#include <Rmath.h>
#include <RcppArmadillo.h>


// =============================================================================
// BaseParameterPrior — abstract base for real-valued parameter priors
//
// Follows the same polymorphic pattern as BaseEdgePrior. Each subclass stores
// its own hyperparameters and provides logp/grad evaluated at a point x.
//
// Used for interaction parameters, threshold parameters, and continuous means.
// =============================================================================
class BaseParameterPrior {
public:
    virtual ~BaseParameterPrior() = default;

    /** Log-density log p(x) up to an additive constant. */
    virtual double logp(double x) const = 0;

    /**
     * Log-density with an additional multiplicative scale factor.
     *
     * For priors with a scale parameter, this evaluates the prior at x
     * with the scale multiplied by scale_factor. Used by OMRF/mixed MRF
     * where the prior scale is adjusted per variable pair based on score
     * range (pairwise_scaling_factors).
     *
     * Default: ignores scale_factor (delegates to logp(x)).
     */
    virtual double logp(double x, double scale_factor) const {
        (void)scale_factor;
        return logp(x);
    }

    /** Gradient d/dx log p(x). */
    virtual double grad(double x) const = 0;

    /**
     * Gradient with an additional multiplicative scale factor.
     * Default: ignores scale_factor (delegates to grad(x)).
     */
    virtual double grad(double x, double scale_factor) const {
        (void)scale_factor;
        return grad(x);
    }

    /**
     * Edge-aware log-density. Carries edge coordinates `(i, j)` for priors
     * whose scale is per-edge (e.g., the Graphical G-prior, where the slab
     * variance is `g · V_ij` with `V_ij` computed from null-MLE Fisher info).
     * Default: ignores (i, j) and delegates to `logp(x)`.
     *
     * Call sites in `GGMModel` and `GGMGradientEngine` pass `(i, j)` when
     * known so an edgewise prior can resolve its per-edge variance; priors
     * with a single shared scale (Cauchy, Normal, …) are unaffected.
     */
    virtual double logp(double x, int i, int j) const {
        (void)i; (void)j;
        return logp(x);
    }

    /** Edge-aware gradient (mirror of edge-aware logp). */
    virtual double grad(double x, int i, int j) const {
        (void)i; (void)j;
        return grad(x);
    }

    /** Deep copy for parallel chains. */
    virtual std::unique_ptr<BaseParameterPrior> clone() const = 0;
};


// =============================================================================
// CauchyPrior — Cauchy(0, scale)
// =============================================================================
class CauchyPrior final : public BaseParameterPrior {
public:
    explicit CauchyPrior(double scale) : scale_(scale) {}

    double logp(double x) const override {
        return R::dcauchy(x, 0.0, scale_, true);
    }

    double logp(double x, double scale_factor) const override {
        return R::dcauchy(x, 0.0, scale_ * scale_factor, true);
    }

    double grad(double x) const override {
        double s2 = scale_ * scale_;
        return -2.0 * x / (s2 + x * x);
    }

    double grad(double x, double scale_factor) const override {
        double s = scale_ * scale_factor;
        double s2 = s * s;
        return -2.0 * x / (s2 + x * x);
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<CauchyPrior>(*this);
    }

private:
    double scale_;
};


// =============================================================================
// NormalPrior — Normal(0, scale)
// =============================================================================
class NormalPrior final : public BaseParameterPrior {
public:
    explicit NormalPrior(double scale) : scale_(scale) {}

    double logp(double x) const override {
        return R::dnorm(x, 0.0, scale_, true);
    }

    double logp(double x, double scale_factor) const override {
        return R::dnorm(x, 0.0, scale_ * scale_factor, true);
    }

    double grad(double x) const override {
        double s2 = scale_ * scale_;
        return -x / s2;
    }

    double grad(double x, double scale_factor) const override {
        double s = scale_ * scale_factor;
        double s2 = s * s;
        return -x / s2;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<NormalPrior>(*this);
    }

private:
    double scale_;
};


// =============================================================================
// BetaPrimePrior — logit-Beta(alpha, beta) prior
//
// If sigma(x) ~ Beta(alpha, beta), then x = logit(Y) where Y ~ Beta(a, b).
// log p(x) = alpha * x - (alpha + beta) * log(1 + exp(x)) + const
// =============================================================================
class BetaPrimePrior final : public BaseParameterPrior {
public:
    BetaPrimePrior(double alpha, double beta)
        : alpha_(alpha), beta_(beta) {}

    double logp(double x) const override {
        return x * alpha_ - std::log1p(std::exp(x)) * (alpha_ + beta_);
    }

    double grad(double x) const override {
        // alpha - (alpha + beta) * sigmoid(x)
        double p = 1.0 / (1.0 + std::exp(-x));
        return alpha_ - (alpha_ + beta_) * p;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<BetaPrimePrior>(*this);
    }

private:
    double alpha_;
    double beta_;
};


// =============================================================================
// GammaScalePrior — Gamma(shape, rate) prior for positive parameters
//
// Used for precision matrix diagonal elements.
// log p(x) = (shape - 1) * log(x) - rate * x + const
// =============================================================================
class GammaScalePrior final : public BaseParameterPrior {
public:
    GammaScalePrior(double shape, double rate)
        : shape_(shape), rate_(rate) {}

    double logp(double x) const override {
        return R::dgamma(x, shape_, 1.0 / rate_, true);
    }

    double grad(double x) const override {
        // d/dx log Gamma(x; shape, rate) = (shape - 1) / x - rate
        return (shape_ - 1.0) / x - rate_;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<GammaScalePrior>(*this);
    }

private:
    double shape_;
    double rate_;
};


// =============================================================================
// GraphicalGPrior — edgewise Normal slab for the Graphical G-prior on GGM
// edge parameters.
//
// Slab:  ω_ij | γ_ij = 1, g  ~  N(0,  g · V_ij)
//
// with V_ij = 1 / (4 n S̄_ii S̄_jj) (Fisher-info-based, computed at fit
// time from the data sufficient statistic). The shared scale parameter
// g = t² is a state variable of the chain; the diagonal Gamma rate is
// tied to 1/t via the companion `GraphicalGDiag` so the scale-matching
// identity Z(Γ, g; δ) = g^{q δ / 2} · Z̃(Γ) holds. See
//   ~/SV/Graphical G-Prior/notes/gg-prior-bgms-implementation-plan.md
// and IDEA.md §5.1 for the construction.
//
// State pointers (V_ij table, t value) are non-owning references back to
// the GGMModel that owns this prior. After cloning the model (parallel
// chains), `bind_to(...)` MUST be called on the clone's prior to rebind
// pointers to the clone's V_ij_table / t_ state. Failing to rebind leaves
// the prior reading the original model's state — a soft aliasing bug.
// =============================================================================
class GraphicalGPrior final : public BaseParameterPrior {
public:
    GraphicalGPrior(const arma::mat* V_ij, const double* t_ptr)
        : V_ij_(V_ij), t_(t_ptr) {}

    /** Rebind state pointers (used by GGMModel copy constructor). */
    void bind_to(const arma::mat* V_ij, const double* t_ptr) {
        V_ij_ = V_ij;
        t_    = t_ptr;
    }

    /** logp at zero-information value (no edge index): degenerate fallback,
     *  used when an edge-blind caller invokes the prior. Treats the slab as
     *  if V_ij ≡ 1, surfacing a defensive default rather than silently using
     *  an arbitrary edge's variance. */
    double logp(double x) const override {
        const double s = *t_;
        return R::dnorm(x, 0.0, s, true);
    }

    double logp(double x, int i, int j) const override {
        const double v = (*V_ij_)(i, j);
        const double s = *t_ * std::sqrt(v);
        return R::dnorm(x, 0.0, s, true);
    }

    double grad(double x) const override {
        const double s2 = (*t_) * (*t_);
        return -x / s2;
    }

    double grad(double x, int i, int j) const override {
        const double v = (*V_ij_)(i, j);
        const double s2 = (*t_) * (*t_) * v;
        return -x / s2;
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        // Pointers carried verbatim; caller must rebind via `bind_to`
        // after cloning the owning model.
        return std::make_unique<GraphicalGPrior>(*this);
    }

private:
    const arma::mat* V_ij_;
    const double*    t_;
};


// =============================================================================
// GraphicalGDiag — Gamma diagonal prior for the Graphical G-prior with rate
// tied to 1 / √g = 1 / t. Shape is fixed at 1 (Exponential), matching the
// scale-matching identity that makes the GGM joint Z(Γ, g) = g^{q δ / 2} ·
// Z̃(Γ) factor cleanly.
//
// Density (rate parameterisation):  p(x; 1, 1/t) = (1/t) · exp(-x / t)
//
// As with GraphicalGPrior, `t_` is a non-owning reference into the owning
// GGMModel and must be rebound on clone.
// =============================================================================
class GraphicalGDiag final : public BaseParameterPrior {
public:
    explicit GraphicalGDiag(const double* t_ptr) : t_(t_ptr) {}

    void bind_to(const double* t_ptr) { t_ = t_ptr; }

    double logp(double x) const override {
        const double t = *t_;
        return R::dgamma(x, 1.0, t, true);   // dgamma uses scale = 1/rate = t
    }

    double grad(double x) const override {
        (void)x;
        return -1.0 / (*t_);
    }

    std::unique_ptr<BaseParameterPrior> clone() const override {
        return std::make_unique<GraphicalGDiag>(*this);
    }

private:
    const double* t_;
};


// =============================================================================
// Factory functions
// =============================================================================

/**
 * Create a parameter prior from a type string and hyperparameters.
 *
 * @param type    One of "cauchy", "normal", "beta-prime"
 * @param scale   Scale for Cauchy/Normal (ignored for beta-prime)
 * @param alpha   Alpha for beta-prime (ignored for Cauchy/Normal)
 * @param beta    Beta for beta-prime (ignored for Cauchy/Normal)
 */
inline std::unique_ptr<BaseParameterPrior> create_parameter_prior(
    const std::string& type,
    double scale = 1.0,
    double alpha = 0.5,
    double beta = 0.5
) {
    if (type == "cauchy") {
        return std::make_unique<CauchyPrior>(scale);
    } else if (type == "normal") {
        return std::make_unique<NormalPrior>(scale);
    } else if (type == "beta-prime") {
        return std::make_unique<BetaPrimePrior>(alpha, beta);
    }
    Rf_error("Unknown parameter prior type: '%s'", type.c_str());
    return nullptr; // unreachable
}


/**
 * Create a scale prior from a type string and hyperparameters.
 *
 * @param type   One of "gamma", "exponential"
 * @param shape  Shape for Gamma (ignored for exponential, set to 1)
 * @param rate   Rate for Gamma/Exponential
 */
inline std::unique_ptr<BaseParameterPrior> create_scale_prior(
    const std::string& type,
    double shape = 1.0,
    double rate = 1.0
) {
    if (type == "gamma" || type == "exponential") {
        double s = (type == "exponential") ? 1.0 : shape;
        return std::make_unique<GammaScalePrior>(s, rate);
    }
    Rf_error("Unknown scale prior type: '%s'", type.c_str());
    return nullptr; // unreachable
}
