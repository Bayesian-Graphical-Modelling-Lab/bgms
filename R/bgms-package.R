#' bgms: Bayesian Analysis of Networks of Binary and/or Ordinal Variables
#'
#' @description
#' The \code{R} package \strong{bgms} provides tools for Bayesian analysis of
#' the ordinal Markov random field (MRF), a graphical model describing networks
#' of binary and/or ordinal variables \insertCite{MarsmanVandenBerghHaslbeck_2025}{bgms}.
#' The likelihood is approximated via a pseudolikelihood, and Markov chain Monte
#' Carlo (MCMC) methods are used to sample from the corresponding pseudoposterior
#' distribution of model parameters.
#'
#' The main entry points are:
#' \itemize{
#'   \item \strong{bgm}: estimation in a one-sample design.
#'   \item \strong{bgmCompare}: estimation and group comparison in an
#'         independent-sample design.
#' }
#'
#' Both functions support Bayesian effect selection with spike-and-slab priors.
#' \itemize{
#'   \item In one-sample designs, \code{bgm} models the presence or absence of
#'   edges between variables. Posterior inclusion probabilities quantify the
#'   plausibility of each edge and can be converted into Bayes factors for
#'   conditional independence tests.
#'
#'   \item \code{bgm} can also model communities (clusters) of variables. The
#'   posterior distribution of the number of clusters provides evidence for or
#'   against clustering \insertCite{SekulovskiEtAl_2025}{bgms}.
#'
#'   \item In independent-sample designs, \code{bgmCompare} estimates group
#'   differences in edge weights and category thresholds. Posterior inclusion
#'   probabilities quantify the evidence for differences and can be converted
#'   into Bayes factors for parameter equivalence tests
#'   \insertCite{MarsmanWaldorpSekulovskiHaslbeck_2024}{bgms}.
#' }
#'
#' @section Tools:
#' The package also provides:
#' \enumerate{
#'   \item Simulation of response data from MRFs with a Gibbs sampler
#'         (\code{\link{simulate_mrf}}).
#'   \item Posterior estimation and edge selection in one-sample designs
#'         (\code{\link{bgm}}).
#'   \item Posterior estimation and group-difference selection in
#'         independent-sample designs (\code{\link{bgmCompare}}).
#' }
#'
#' @section Vignettes:
#' For tutorials and worked examples, see:
#' \itemize{
#'   \item \code{vignette("intro", package = "bgms")} — Getting started.
#'   \item \code{vignette("comparison", package = "bgms")} — Model comparison.
#'   \item \code{vignette("diagnostics", package = "bgms")} — Diagnostics and
#'         spike-and-slab summaries.
#' }
#'
#' @docType package
#' @keywords internal
#' @useDynLib bgms, .registration=TRUE
#' @references
#' \insertAllCited{}
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
