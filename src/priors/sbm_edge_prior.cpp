#include <RcppArmadillo.h>
#include "rng/rng_utils.h"
#include "math/explog_switch.h"

// ----------------------------------------------------------------------------|
// The c++ code below is based on the R code accompanying the paper:
//  Geng, J., Bhattacharya, A., & Pati, D. (2019). Probabilistic Community
//  Detection With Unknown Number of Communities, Journal of the American
//  Statistical Association, 114:526, 893-905, DOI:10.1080/01621459.2018.1458618
// ----------------------------------------------------------------------------|

// ----------------------------------------------------------------------------|
// A c++ version of table
// ----------------------------------------------------------------------------|
arma::uvec table_cpp(arma::uvec x) {
  arma::uword n = x.n_elem;
  arma::uword m = arma::max(x);
  arma::uvec counts(m+1,arma::fill::zeros);

  for (arma::uword i = 0; i < n; i++) {
    counts(x(i))++;
  }

  return counts;
}


// ----------------------------------------------------------------------------|
// Add a row and column to a matrix (and fill with beta variables)
// Modified to support separate within/between cluster hyperparameters
// ----------------------------------------------------------------------------|
arma::mat add_row_col_block_prob_matrix(arma::mat X,
                                        double beta_alpha,
                                        double beta_beta,
                                        SafeRNG& rng,
                                        double beta_bernoulli_alpha_between,
                                        double beta_bernoulli_beta_between) {
  arma::uword dim = X.n_rows;
  arma::mat Y(dim+1,dim+1,arma::fill::zeros);

  for(arma::uword r = 0; r < dim; r++) {
    for(arma::uword c = 0; c < dim; c++) {
      Y(r, c) = X(r, c);
    }
  }

  // Add new row and column for the new cluster
  for(arma::uword i = 0; i < dim; i++) {
    // Between-cluster edge probabilities (new cluster to existing clusters)
    Y(dim, i) = rbeta(rng, beta_bernoulli_alpha_between, beta_bernoulli_beta_between);
    Y(i, dim) = Y(dim, i);
  }

  // Within-cluster edge probability (diagonal element for new cluster)
  Y(dim, dim) = rbeta(rng, beta_alpha, beta_beta);

  return Y;
}



/**
 * Function: log_likelihood_mfm_sbm
 *
 * Computes the log-likelihood contribution for a single node under the
 * Mixture of Finite Mixtures Stochastic Block Model (MFM-SBM). Evaluates
 * the probability of observed edges between the node and all other nodes
 * given their cluster assignments and cluster connection probabilities.
 *
 * Inputs:
 *  - cluster_assign: Vector of cluster assignments for all nodes.
 *  - cluster_probs: Matrix of edge probabilities between clusters.
 *  - indicator: Symmetric binary matrix of edge indicators (1 = edge present).
 *  - node: Index of the node whose contribution is computed.
 *  - no_variables: Total number of nodes in the network.
 *
 * Returns:
 *  - Log-likelihood contribution for the specified node.
 */
double log_likelihood_mfm_sbm(arma::uvec cluster_assign,
                              arma::mat cluster_probs,
                              arma::umat indicator,
                              arma::uword node,
                              arma::uword no_variables) {
  double output = 0;

  for(arma::uword j = 0; j < no_variables; j++) {
    if(j != node) {
      if(j < node) {
        output += indicator(j, node) *
          MY_LOG(cluster_probs(cluster_assign(j), cluster_assign(node)));
        output += (1 - indicator(j, node)) *
          MY_LOG(1 - cluster_probs(cluster_assign(j), cluster_assign(node)));
      } else {
        output += indicator(node, j) *
          MY_LOG(cluster_probs(cluster_assign(node), cluster_assign(j)));
        output += (1 - indicator(node, j)) *
          MY_LOG(1 - cluster_probs(cluster_assign(node), cluster_assign(j)));
      }
    }
  }

  return output;
}

/**
 * Function: log_marginal_mfm_sbm
 *
 * Computes the log-marginal likelihood contribution for a single node under
 * the MFM-SBM after integrating out cluster connection probabilities. Uses
 * Beta-Bernoulli conjugacy with separate hyperparameters for within-cluster
 * and between-cluster edges.
 *
 * Inputs:
 *  - cluster_assign: Vector of cluster assignments for all nodes.
 *  - indicator: Symmetric binary matrix of edge indicators (1 = edge present).
 *  - node: Index of the node whose contribution is computed.
 *  - no_variables: Total number of nodes in the network.
 *  - beta_bernoulli_alpha: Alpha hyperparameter for within-cluster edges.
 *  - beta_bernoulli_beta: Beta hyperparameter for within-cluster edges.
 *  - beta_bernoulli_alpha_between: Alpha hyperparameter for between-cluster edges.
 *  - beta_bernoulli_beta_between: Beta hyperparameter for between-cluster edges.
 *  - sbm_singleton_boost: If true, add a balancing term for new singleton clusters
 *                         to compensate for the asymmetry in within/between priors.
 *
 * Returns:
 *  - Log-marginal likelihood contribution for the specified node.
 */
double log_marginal_mfm_sbm(arma::uvec cluster_assign,
                            arma::umat indicator,
                            arma::uword node,
                            arma::uword no_variables,
                            double beta_bernoulli_alpha,
                            double beta_bernoulli_beta,
                            double beta_bernoulli_alpha_between,
                            double beta_bernoulli_beta_between,
                            bool sbm_singleton_boost) {

  arma::uvec indices = arma::regspace<arma::uvec>(0, no_variables-1); // vector of variables indices [0, 1, ..., no_variables-1]
  arma::uvec select_variables = indices(arma::find(indices != node)); // vector of variables indices excluding 'node'
  arma::uvec cluster_assign_wo_node = cluster_assign(select_variables); // vector of cluster labels for all variables but excluding 'node'
  arma::uvec indicator_node = indicator.col(node); // column of indicator matrix corresponding to 'node'
  arma::vec gamma_node = arma::conv_to<arma::vec>::from(indicator_node(select_variables)); // selecting only indicators between 'node' and the remaining variables (thus excluding indicator of node with itself -- that is indicator[node,node])
  arma::uvec table_cluster = table_cpp(cluster_assign_wo_node); // frequency table of clusters excluding node

  // Get the cluster assignment of the current node
  arma::uword node_cluster = cluster_assign(node);
  
  // Check if node is in a NEW singleton cluster (not present in other nodes' assignments)
  bool is_new_singleton = (node_cluster >= table_cluster.n_elem);

  double output = 0;
  
  // If sbm_singleton_boost is enabled and this is a new singleton cluster,
  // add a balancing term to compensate for the asymmetry where only between-cluster
  // priors are used. This makes the within-cluster prior influence cluster creation.
  if (sbm_singleton_boost && is_new_singleton) {
    double within_prior_mean = beta_bernoulli_alpha / (beta_bernoulli_alpha + beta_bernoulli_beta);
    double between_prior_mean = beta_bernoulli_alpha_between / (beta_bernoulli_alpha_between + beta_bernoulli_beta_between);
    
    // Only add boost if within-prior favors more connectivity than between-prior
    if (within_prior_mean > between_prior_mean) {
      // Count total edges from this node
      double total_edges = arma::accu(gamma_node);
      double total_pairs = static_cast<double>(no_variables - 1);
      double edge_rate = total_edges / total_pairs;
      
      // Compute the log-likelihood ratio of observing these edges under
      // within-cluster prior vs between-cluster prior.
      // This represents how much more favorable the within-cluster prior
      // would be for this node's connectivity pattern.
      double log_ratio_edges = 0.0;
      if (edge_rate > 0 && edge_rate < 1) {
        log_ratio_edges = total_pairs * (
          edge_rate * std::log(within_prior_mean / between_prior_mean) +
          (1.0 - edge_rate) * std::log((1.0 - within_prior_mean) / (1.0 - between_prior_mean))
        );
      }
      
      // Scale down to avoid overcompensating - this is a heuristic adjustment
      // to partially offset the between-cluster penalty without fully replacing it.
      // The 0.5 factor represents uncertainty about future cluster growth.
      output += 0.5 * std::max(0.0, log_ratio_edges);
    }
  }
  
  for(arma::uword i = 0; i < table_cluster.n_elem; i++){
    if(table_cluster(i) > 0){ // if the cluster is empty -- table_cluster(i) = 0 == then it is the previous cluster of 'node' where 'node' was the only member - a singleton, thus skip)
      arma::uvec which_variables_cluster_i = arma::find(cluster_assign_wo_node == i); // which variables belong to cluster i
      double sumG = arma::accu(gamma_node(which_variables_cluster_i)); // sum the indicator variables between node and those variables
      double sumN = static_cast<double>(which_variables_cluster_i.n_elem); // take the size of the group as maximum number of relations

      // Determine if this is within-cluster or between-cluster
      bool is_within_cluster = (i == node_cluster);

      double alpha = is_within_cluster ? beta_bernoulli_alpha : beta_bernoulli_alpha_between;
      double beta = is_within_cluster ? beta_bernoulli_beta : beta_bernoulli_beta_between;

      output += R::lbeta(sumG + alpha, sumN - sumG + beta) - R::lbeta(alpha, beta); // calculate log-density for cluster i and sum it to the marginal log-likelihood
    }
  }
  return output;
}

// ----------------------------------------------------------------------------|
// Helper function to update sumG in sample_block_probs_mfm_sbm()
// ----------------------------------------------------------------------------|
inline void update_sumG(double &sumG,
                        const arma::uvec &cluster_assign,
                        const arma::umat &indicator,
                        arma::uword r,
                        arma::uword s,
                        arma::uword no_variables) {
  for(arma::uword node1 = 0; node1 < no_variables - 1; node1++) {
    if(cluster_assign(node1) == r) {
      for(arma::uword node2 = node1 + 1; node2 < no_variables; node2++) {
        if(cluster_assign(node2) == s) {
          sumG += static_cast<double>(indicator(node1, node2));
        }
      }
    }
  }
}

// ----------------------------------------------------------------------------|
// Sample the cluster assignment in sample_block_allocations_mfm_sbm()
// ----------------------------------------------------------------------------|
arma::uword sample_cluster(arma::vec cluster_prob,
                           SafeRNG& rng) {
  arma::vec cum_prob = arma::cumsum(cluster_prob);
  double u = runif(rng) * arma::max(cum_prob);

  for (arma::uword i = 0; i < cum_prob.n_elem; i++) {
    if (u <= cum_prob(i)) {
      return i;
    }
  }
  return cum_prob.n_elem;
}

// ----------------------------------------------------------------------------|
// Sample the block allocations for the MFM - SBM
// Modified to support separate within/between cluster hyperparameters
// ----------------------------------------------------------------------------|
arma::uvec block_allocations_mfm_sbm(arma::uvec cluster_assign,
                                     arma::uword no_variables,
                                     arma::vec log_Vn,
                                     arma::mat block_probs,
                                     arma::umat indicator,
                                     arma::uword dirichlet_alpha,
                                     double beta_bernoulli_alpha,
                                     double beta_bernoulli_beta,
                                     double beta_bernoulli_alpha_between,
                                     double beta_bernoulli_beta_between,
                                     bool sbm_singleton_boost,
                                     SafeRNG& rng) {

  arma::uword old;
  arma::uword cluster;
  arma::uword no_clusters;
  double prob;
  double loglike;
  double logmarg;

  // Generate a randomized order using Rcpp's sample function
  arma::uvec indices = arma_randperm(rng, no_variables);

  for (arma::uword idx = 0; idx < no_variables; idx++) {
    arma::uword node = indices(idx);
    old = cluster_assign(node);

    arma::uvec cluster_size = table_cpp(cluster_assign);
    no_clusters = cluster_size.n_elem;

    if (cluster_size(old) == 1) {
      // Singleton cluster.

      // Cluster sizes without node
      arma::uvec cluster_size_node = cluster_size;

      // Compute probabilities for sampling process
      arma::vec cluster_prob(no_clusters + 1);
      for (arma::uword c = 0; c <= no_clusters; c++) {
        arma::uvec cluster_assign_tmp = cluster_assign;
        cluster_assign_tmp(node) = c;

        if (c < no_clusters) {
          if(c != old){
            loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                             block_probs,
                                             indicator,
                                             node,
                                             no_variables);

            prob = (static_cast<double>(dirichlet_alpha) + static_cast<double>(cluster_size_node(c))) *
              MY_EXP(loglike);
          }
          else{ // if old group, the probability is set to 0.0
            prob = 0.0;
          }

        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         indicator,
                                         node,
                                         no_variables,
                                         beta_bernoulli_alpha,
                                         beta_bernoulli_beta,
                                         beta_bernoulli_alpha_between,
                                         beta_bernoulli_beta_between,
                                         sbm_singleton_boost);

          prob = static_cast<double>(dirichlet_alpha) *
            MY_EXP(logmarg) *
            MY_EXP(log_Vn(no_clusters - 1) - log_Vn(no_clusters - 2));
        }

        cluster_prob(c) = prob;
      }

      //Choose the cluster number for node
      cluster = sample_cluster(cluster_prob, rng);

      //if the sampled cluster is the new added cluster or the old one
      if (cluster == no_clusters) {
        cluster_assign(node) = old; // new cluster takes the place of the older singleton but doesn't update probabilities, they are kept the same as the old ones
      } else { // otherwise remove old (singleton) empty cluster, redefine cluster_assign and block_probs
        cluster_assign(node) = cluster;
        for (arma::uword i = 0; i < no_variables; i++) {
          if (cluster_assign(i) > old) {
            cluster_assign(i) -= 1;
          }
        }
        // removing row and col index 'old' from block_probs
        block_probs.shed_row(old);
        block_probs.shed_col(old);
      }
    } else {
      // Cluster sizes without node
      arma::uvec cluster_size_node = cluster_size;
      cluster_size_node(old) -= 1;

      // Compute probabilities for sampling process
      arma::vec cluster_prob(no_clusters + 1);
      for (arma::uword c = 0; c <= no_clusters; c++) {
        arma::uvec cluster_assign_tmp = cluster_assign;
        cluster_assign_tmp(node) = c;
        if (c < no_clusters) {
          loglike = log_likelihood_mfm_sbm(cluster_assign_tmp,
                                           block_probs,
                                           indicator,
                                           node,
                                           no_variables);

          prob = (static_cast<double>(dirichlet_alpha) + static_cast<double>(cluster_size_node(c))) *
            MY_EXP(loglike);
        } else {
          logmarg = log_marginal_mfm_sbm(cluster_assign_tmp,
                                         indicator,
                                         node,
                                         no_variables,
                                         beta_bernoulli_alpha,
                                         beta_bernoulli_beta,
                                         beta_bernoulli_alpha_between,
                                         beta_bernoulli_beta_between,
                                         sbm_singleton_boost);

          prob = static_cast<double>(dirichlet_alpha) *
            MY_EXP(logmarg) *
            MY_EXP(log_Vn(no_clusters) - log_Vn(no_clusters-1));
        }

        cluster_prob(c) = prob;
      }


      //Choose the cluster number for node
      cluster = sample_cluster(cluster_prob, rng);

      cluster_assign(node) = cluster;

      if (cluster == no_clusters) {
        block_probs = add_row_col_block_prob_matrix(block_probs,
                                                    beta_bernoulli_alpha,
                                                    beta_bernoulli_beta,
                                                    rng,
                                                    beta_bernoulli_alpha_between,
                                                    beta_bernoulli_beta_between);
      }
    }
  }
  return cluster_assign;

}

// ----------------------------------------------------------------------------|
// Sample the block parameters for the MFM - SBM
// Modified to support separate within/between cluster hyperparameters
// ----------------------------------------------------------------------------|
arma::mat block_probs_mfm_sbm(arma::uvec cluster_assign,
                              arma::umat indicator,
                              arma::uword no_variables,
                              double beta_bernoulli_alpha,
                              double beta_bernoulli_beta,
                              double beta_bernoulli_alpha_between,
                              double beta_bernoulli_beta_between,
                              SafeRNG& rng) {

  arma::uvec cluster_size = table_cpp(cluster_assign);
  arma::uword no_clusters = cluster_size.n_elem;

  arma::mat block_probs(no_clusters, no_clusters);
  arma::mat theta(no_variables, no_variables);

  double sumG;
  double size;

  for(arma::uword r = 0; r < no_clusters; r++) {
    for(arma::uword s = r; s < no_clusters; s++) {
      sumG = 0;

      if(r == s) {
        // Within-cluster: always use main parameters
        update_sumG(sumG, cluster_assign, indicator, r, r, no_variables);
        size = static_cast<double>(cluster_size(r)) * (static_cast<double>(cluster_size(r)) - 1) / 2;
        block_probs(r, s) = rbeta(rng,
                    sumG + beta_bernoulli_alpha,
                    size - sumG + beta_bernoulli_beta);
      } else {
        // Between-cluster: use between parameters
        update_sumG(sumG, cluster_assign, indicator, r, s, no_variables);
        update_sumG(sumG, cluster_assign, indicator, s, r, no_variables);
        size = static_cast<double>(cluster_size(s)) * static_cast<double>(cluster_size(r));

        block_probs(r, s) = rbeta(rng, sumG + beta_bernoulli_alpha_between, size - sumG + beta_bernoulli_beta_between);
      }
      block_probs(s, r) = block_probs(r, s);
    }
  }

  return block_probs;
}