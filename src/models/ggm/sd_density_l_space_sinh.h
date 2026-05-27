#pragma once

#include <RcppArmadillo.h>

// Sinh-substitution + midpoint-rule quadrature for the L-space SD primitive.
//
// Computes the normaliser of the conditional density
//
//   I(A, B, s, alpha) = integral over R of (s + phi^2)^(alpha-1) *
//                       exp(-A phi^2 + B phi) dphi
//
// via the substitution
//
//   phi = sqrt(s) * sinh(t)
//
// which transforms the integral to
//
//   I = s^(alpha - 1/2) * integral over R of g(t) dt
//
//   g(t) = cosh(t)^(2 alpha - 1) *
//          exp(-A s sinh^2(t) + B sqrt(s) sinh(t))
//
// Structural advantage. In the phi coordinate the branch points of
// (s + phi^2)^(alpha-1) sit at phi = +- i sqrt(s) -- distance sqrt(s)
// from the real axis, which collapses for small s. The sinh substitution
// maps these to t = +- i pi/2 (fixed location). For alpha > 1/2 the
// transformed integrand has zeros there; for alpha <= 1/2 it has an
// integrable singularity. In every case the distance from the real axis
// is bounded below by pi/2, so the quadrature convergence rate becomes
// uniform across the (alpha, s) plane.
//
// Quadrature. Midpoint rule on a truncated interval [T_left, T_right]
// centred at the Gaussian-envelope mode in t (t* = arcsinh(B/(2 A
// sqrt(s)))). The interval covers +- n_sigma standard deviations of the
// envelope in the sinh coordinate (Gaussian-in-sinh variance = 1/(2 A s)).
//
// Convergence. For midpoint rule on a truncated interval [T_left, T_right]
// with N nodes, applied to an analytic integrand on R with the nearest
// complex singularity at distance b from the real axis, the error decays
// as exp(- pi * b * N / (T_right - T_left)) up to polynomial prefactors.
// Here b >= pi/2 so the convergence rate is robust.
//
// Status codes:
//   0  ok.
//   1  A <= 0 (PD revert condition).
//   2  s_jj <= 0 (PD revert condition).

namespace ggm_sd {

struct LSDSinhResult {
    double log_density;  ///< log pi(x_eval | rest, Y)
    double log_Z;        ///< log of the normaliser
    int    status;
};

LSDSinhResult density_at_l_ji_sinh(double x_eval,
                                    double A, double B,
                                    double s_jj,
                                    double alpha,
                                    int num_nodes = 128);

}  // namespace ggm_sd
