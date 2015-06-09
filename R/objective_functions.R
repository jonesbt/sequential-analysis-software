#' Computes the coefficient of variance for a beta distribution with shape
#' parameters alpha and beta - alpha. This function is private to this file.
cv.beta <- function(alpha, beta) {
    if(beta < alpha)
        stop('beta must be greater than alpha for a real answer.')
    if(alpha == 0 || beta <= 1)
        return(NA)
    return(sqrt((beta - alpha) / 
                (alpha * (beta - 1))))
}

#' Computes 1 - f(delta, alpha_i, alpha_0), where f(delta, alpha_i, alpha_0)
#' is the CDF of the beta distribution with shape parameters alpha_i and alpha_0
#' evaluated at delta. This function is private to this file.
compute_prob = function(alpha_0, alpha_i, delta) {
    return(1 - pbeta(delta, alpha_i, alpha_0))
}

#' Compute the objective function that is described in the text.
#'
#' \param alphas An m-by-k matrix for a system with m origins and k
#' destinations. alphas[i,j] is x[i,j] as described in the manuscript text.
#' 
#' The following arguments should be in obj_fn_args.
#' \param delta The threshold for "irrelevant" connections.
#' \param epsilon The threshold for stopping.
#' \param pi The probability that a connection must be "irrelevant" for it to
#' be ignored.
#' \return The value of the objective function.
obj_fn_probabilities = function(alphas, obj_fn_args) {
    delta = obj_fn_args$delta
    epsilon = obj_fn_args$epsilon
    pi = obj_fn_args$pi
    # Compute the C.V. for each p_{ij}.
    cvs = matrix(NA, nrow(alphas), ncol(alphas))
    for(i in seq(nrow(alphas)))
        for(j in seq(ncol(alphas)))
            if((alphas[i,j] / sum(alphas[i,])) > delta)
                cvs[i,j] = cv.beta(alphas[i,j], sum(alphas[i,]))
    # Compute the probability that each p_{ij} > delta.
    probs = matrix(NA, nrow(alphas), ncol(alphas))
    for(i in seq(nrow(alphas)))
        for(j in seq(ncol(alphas)))
            if((alphas[i,j] / sum(alphas[i,])) < delta)
                probs[i,j] = compute_prob(delta, sum(alphas[i,]), alphas[i,j])
    probs[is.na(probs)] = 1
    cvs[cvs < epsilon] = 0 # Condition 1
    cvs[probs < pi] = 0 # Condition 2
    # Compute the cost for each site.
    site_costs = rep(NA, nrow(alphas))
    for(i in seq(nrow(alphas))) {
        if(sum(is.na(cvs[i,])) == length(cvs[i,])) {
            max_cv = 0
        } else {
            max_cv = max(cvs[i,], na.rm=TRUE)
        }
        site_costs[i] = max_cv
    }
    # Now return the result.
    return(max(site_costs))    
}
