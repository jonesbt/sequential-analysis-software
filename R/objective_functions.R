#' Compute the mean of a beta distribution with parameters alpha and beta. This
#' function is private to this file.
#' @param alpha The first shape parameter of the beta distribution.
#' @param beta The second shape parameter of the beta distribution.
#' @return The expected value for the beta distribution.
mu.beta <- function(alpha, beta) {
    return(alpha / (alpha + beta))
}

#' Compute the standard deviation of a beta distribution with parameters alpha
#' and beta. This function is private to this file.
#' @param alpha The first shape parameter of the beta distribution.
#' @param beta The second shape parameter of the beta distribution.
#' @return The standard deviation for the beta distribution.
sd.beta <- function(alpha, beta) {
    return(sqrt((alpha * beta) /
                ((alpha + beta)^2 * (alpha + beta + 1))))
}

#' Computes the coefficient of variance for a beta distribution with
#' parameters alpha and beta. This function is private to this file.
#' @param alpha The first shape parameter of the beta distribution.
#' @param beta The second shape parameter of the beta distribution.
#' @return The standard deviation for the beta distribution.
cv.beta <- function(alpha, beta) {
    return(sd.beta(alpha, beta) / mu.beta(alpha, beta))
}

#' Computes 1 - f(delta, alpha, beta), where f(delta, alpha, beta)
#' is the CDF of the beta distribution with shape parameters alpha and beta
#' evaluated at delta. This function is private to this file.
compute_prob = function(delta, alpha, beta) {
    return(1 - pbeta(delta, alpha, beta))
}

#' Compute the objective function that is described in the text.
#'
#' @param alphas An m-by-k matrix for a system with m origins and k
#' destinations. alphas[i,j] is x[i,j] as described in the manuscript text.
#' 
#' The following arguments should be in obj_fn_args.
#' @param delta The threshold for "irrelevant" connections.
#' @param epsilon The threshold for stopping.
#' @param pi The probability that a connection must be "irrelevant" for it to
#' be ignored.
#' @return The value of the objective function.
#' @export
obj_fn_probabilities = function(alphas, obj_fn_args, by_site=FALSE) {
    delta = obj_fn_args$delta
    epsilon = obj_fn_args$epsilon
    pi = obj_fn_args$pi
    # Compute the C.V. for each p_{ij}.
    cvs = matrix(NA, nrow(alphas), ncol(alphas))
    for(i in seq(nrow(alphas)))
        for(j in seq(ncol(alphas)))
            if((alphas[i,j] / sum(alphas[i,])) > delta)
                cvs[i,j] = cv.beta(alphas[i,j], sum(alphas[i,]) - alphas[i,j])
    # Compute the probability that each p_{ij} > delta.
    probs = matrix(NA, nrow(alphas), ncol(alphas))
    for(i in seq(nrow(alphas)))
        for(j in seq(ncol(alphas)))
            if((alphas[i,j] / sum(alphas[i,])) < delta)
                probs[i,j] = compute_prob(delta, alphas[i,j],
                                          sum(alphas[i,]) - alphas[i,j])
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
    if(by_site)
        return(site_costs)
    return(max(site_costs))    
}
