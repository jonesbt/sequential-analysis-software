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
#' @param n_sites The number of sites. This should be a perfect square
#' (4, 9, ...).
#' @param topology_params A list with elements p and nei, which are passed
#' through to igraph::watts.strogatz.game().
#' @export
generate_conn_mat = function(n_sites, topology_params=NULL) {
    # Generate the adjacency matrix.
    adj_mat = generate_adj_mat_small_world(n_sites, topology_params)
    # Generate the connectivity matrix.
    return(generate_conn_mat_unif(adj_mat))
}

#' \todo Add test that this is accurately simulating the matrix.
sample_conn_mat = function(conn_mat, origin, n) {
    return(apply(rmultinom(n, 1, conn_mat[origin,]), 2,
                 function(x) which(x == 1)))
}

#' @export
normalize_distribution = function(dist, n) {
    dist = floor(dist)
    i = 1
    while(sum(dist) < n) {
        dist[i] = dist[i] + 1
        i = i + i
        if(i > length(dist))
            i = 1
    }
    return(dist)
}

#' @export
update_alphas = function(alphas, data, start, count) {
    # Compute the ending sample to use.
    end = start + count - 1
    # Update alphas in place.
    for(i in seq(nrow(alphas)))
        if(count[i] > 0)
            for(j in seq(ncol(alphas)))
                alphas[i,j] = alphas[i,j] + sum(data[i,start[i]:end[i]] == j)
    # Return the results including the next starting count.
    return(list(alphas=alphas, start=end + 1))
}
