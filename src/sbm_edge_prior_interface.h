#include <RcppArmadillo.h>
#include "math/explog_switch.h"



arma::vec compute_Vn_mfm_sbm(arma::uword num_variables,
                             double dirichlet_alpha,
                             arma::uword t_max,
                             double lambda);