#' Computes the expected posterior value of the objective function.
#' This function is private to this file.
#' 
#' \param alphas A row vector from the matrix alphas given in the text.
#' \param n The number of new tracers to release.
#' \param obj_fn The objective function.
#' \param obj_fn_args Arguments to pass to the objective function.
#' \return The expected posterior value of the objective function.
compute_exp_cost = function(alphas, n, obj_fn, obj_fn_args) {
    ## Save the original parameters.
    orig_alphas = matrix(alphas, nrow=1)
    n_reps = 2.5e2
    # Compute the C.V. for each p_{ij}.
    costs = numeric(n_reps)
    for(r in seq(n_reps)) {
        ## Restore the original values of alpha.
        alphas = orig_alphas
        ## Generate a probability vector.
        p = rdirichlet(1, alphas)
        ## Generate a random sample for where these particles go.
        dest = rmultinom(1, n, p)
        ## Update the alphas.
        alphas = alphas + t(dest)
        # Now compute the cost on the updated alphas.
        costs[r] = obj_fn(alphas, obj_fn_args)
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
    ## Initially create an empty distribution.
    release_distribution = integer(nrow(alphas))
    ## Create a vector to store the costs and the allocation info.
    release_site = NA
    ## Compute the initial costs.
    prior_costs = obj_fn_probabilities(alphas, by_site=TRUE, obj_fn_args)
    exp_costs = sapply(seq(nrow(alphas)), function(j)
        compute_exp_cost(alphas[j,], block_size, obj_fn, obj_fn_args))
    ## Iterate through the particles.
    for(i in seq(1, n, by=block_size)) {
        ## Update the costs if necessary.
        if(!is.na(release_site)) {
            prior_costs[release_site] = exp_costs[release_site]
            exp_costs[release_site] =
                compute_exp_cost(alphas[release_site,],
                                 release_distribution[release_site] +block_size,
                                 obj_fn, obj_fn_args)
        }
        costs = sapply(seq(length(exp_costs)), function(i) {
            x = apply(cbind(prior_costs, exp_costs), 1, max)
            x[i] = exp_costs[i]
            return(max(x))
        })
        release_site = which(costs == min(costs))
        if(length(release_site) > 1) {
            release_site = which.max(prior_costs)
        }
        ## Release the particle from the site that minimizes the expected cost.
        release_distribution[release_site] =
            release_distribution[release_site] + block_size
    }
    ## Uniformly allocate any remaining particles.
    i = 1
    while(sum(release_distribution) < n) {
        release_distribution[i] = release_distribution[i] + 1
        i = i + 1
        if(i > length(release_distribution))
            i = 1
    }
    ## Return the result.
    return(list(dist=release_distribution, cost=max(costs)))
}
