


# approx_f_xy() : mixture density used to approximate the true density,
#                 covariate-INDEPENDENT weights + variable selection
# input:  
#          theta  : variational parameters
#          prior  : prior-related variables
#          y      : (1 x 1) response
#          x      : (1 x D) covariates
# output: 
#          p_y    : approximate conditional density for y|X 
approx_f_xy = function(prior, theta, y, x) {
    
    
    #if (theta$VERBOSE) {
    #    print("MoE approx -- covariate-INDEPENDENT weights + V.S.")
    #}
    
    p_y = 0
    
    # print(x)
    
    # evaluate each of the gaussians
    for (k in 1:prior$K) {
        tau_k_inv = 1 / theta$tau_k[k]  # precision of k-th gaussian
        beta_k = theta$beta_k[,k]         # coefficient vector of k-th gaussian
        mu_k = c(t(x) %*% beta_k)       # mean component of k-th gaussian
        # print(mu_k)
        p_y = p_y + theta$pi_k * dnorm(y, mu_k, sqrt(tau_k_inv))
    }
    
    # print(p_y)
    
    return(p_y)
} # end of approx_f_xy() function


