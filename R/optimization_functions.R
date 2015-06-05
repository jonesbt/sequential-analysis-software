#' Computes the expected posterior value of the objective function.
#' This function is private to this file.
#' 
#' \param alphas A row vector from the matrix alphas given in the text.
#' \param n The number of new tracers to release.
#' \param obj_fn The objective function.
#' \param obj_fn_args Arguments to pass to the objective function.
#' \return The expected posterior value of the objective function.
compute_exp_cost = function(alphas, n, obj_fn, obj_fn_args) {
    alphas = matrix(alphas, nrow=1)
    orig_alphas = alphas
    n_reps = 1e2
    # Compute the C.V. for each p_{ij}.
    costs = numeric(n_reps)
    for(i in seq(n_reps)) {
        alphas = orig_alphas
        # Generate a random sample for where these particles go.
        dest = rmultinom(1, n, alphas / sum(alphas))
        alphas = alphas + t(dest)
        # Now compute the cost on the updated alphas.
        costs[i] = obj_fn(alphas, obj_fn_args)
    }
    return(mean(costs))
}

#' This function implements a greedy optimization heuristic. At each step, we
#' compute the expected posterior value of the objective function after
#' allocating the particle to each node. We then assign the particle to the
#' site with the lowest of these values.
#'
#' \param alphas The matrix alpha as described in the manuscript text.
#' \param n The number of new tracers to allocate.
#' \param obj_fn The objective function to optimize over.
#' \param obj_fn_args Arguments to pass to the objective function.
#' \param block_size An optional number of tracers to allocate as a single
#' block. The algorithm is approximately linear with respect to the number of
#' blocks, so setting a block size greater than 1 can be useful for reducing
#' runtime with very large number of tracers.
#' \return A list with the release distribution (dist) and the expected
#' posterior cost from using this distribution (cost).
optimization_heuristic = function(alphas, n, 
    obj_fn, obj_fn_args, block_size=1) {
    # Initially create an empty distribution.
    release_distribution = integer(nrow(alphas))
    # Compute the costs for the first iteration.
    costs = sapply(seq(nrow(alphas)), function(j)
        compute_exp_cost(alphas[j,], block_size, obj_fn, obj_fn_args)
    # Create the prior connectivity matrix. This is updated every timestep.
    P = alphas
    for(i in seq(1, n, by=block_size)) {
        # Release the particle from the site with the highest cost.
        release_site = which.max(costs)
        release_distribution[release_site] =
            release_distribution[release_site] + block_size
        # Normalize the probabilities and update the alphas based on the
        # expected cost.
        for(j in seq(nrow(P)))
            P[release_site,] = (alphas[release_site,] ) /
                sum(alphas[release_site,] )
        alphas[release_site, ] = alphas[release_site, ] +
            P[release_site, ]
        costs[release_site] = compute_exp_cost(alphas[release_site,],
             block_size, obj_fn, obj_fn_args)
    }
    i = 1
    while(sum(release_distribution) < n) {
        release_distribution[i] = release_distribution[i] + 1
        i = i + 1
        if(i > length(release_distribution))
            i = 1
    }
    return(list(dist=release_distribution, cost=max(costs)))
}
