# This file contains functions for generating connectivity matrices.

#' Generates a small world topology. This function is private to this file.
generate_adj_mat_small_world = function(n_sites, params) {
  n = sqrt(n_sites)
  if(n != floor(n))
    stop('The number of sites must have an integer square root for small world
      matrices.')
  if(is.null(params$nei))
    stop('Must supply nei parameter.')
  if(is.null(params$p))
    stop('Must supply p parameter.')
  g = watts.strogatz.game(2, n, params$nei, params$p)
  return(get.adjacency(g, type='both', sparse=FALSE) + diag(n_sites))
}

#' Divides the connection strengths among the edges as described in the text.
#' This function is private to this file.
generate_conn_mat_unif = function(adj_mat) {
  n_sites = nrow(adj_mat) 
  for(i in seq(nrow(adj_mat))) {
    prop_success = runif(1, min=0, max=0.1)
    # Now generate n_dest uniform random numbers that sum to prop_succes.
    n_dest = sum(adj_mat[i,])
    vals = rev(sort(c(0, runif(n_dest - 1, 0, prop_success))))
    vals = vals[-length(vals)] - vals[-1]
    # Make those number the edge weights.
    adj_mat[i, adj_mat[i,] > 0] = rdirichlet(1, rep(1, n_dest)) * prop_success
  }
  adj_mat = cbind(adj_mat, 1 - apply(adj_mat, 1, sum))
  return(adj_mat)
}

#' Generates a connectivity matrix as described in the manuscript text.
#'
#' \param n_sites The number of sites. This should be a perfect square
#' (4, 9, ...).
#' \param topology_params A list with elements p and nei, which are passed
#' through to igraph::watts.strogatz.game().
generate_conn_mat = function(n_sites, topology_params=NULL) {
    # Generate the adjacency matrix.
    adj_mat = generate_adj_mat_small_world(n_sites, topology_params)
    # Generate the connectivity matrix.
    return(generate_conn_mat_unif(adj_mat))
}
