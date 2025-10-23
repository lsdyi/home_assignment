RWMsampler <- function(logPostFunc, initVal, nSim, nBurn, Sigma, c, ...) {
    # Run the algorithm for nSim iterations
    # using the multivariate proposal N(theta_previous_draw, c*Sigma)
    # Return the posterior draws after discarding nBurn iterations as burn-in
    draws <- matrix(rep(0, nSim * length(initVal)), ncol = length(initVal))
    for (i in 1:nSim) {
        temp <- as.vector(rmvnorm(1, mean = initVal, sigma = c * Sigma))
        alpha <- min(1, exp(logPostFunc(temp, ...) - logPostFunc(initVal, ...)))
        u <- runif(1)
        if (alpha >= u) {
            # effective draw
            draws[i, ] <- temp
        } else {
            draws[i, ] <- initVal
        }

        # update previous value
        initVal <- draws[i, ]
    }

    return(draws[(nBurn + 1):nSim, ])
}
