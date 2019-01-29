

## approxDensity.R
## generate the approximating density using a mixture of experts, using
## the model parameters from CAVI


# input:  
#          theta  : 
#          prior  : 
#          K      :
#          data   :
# output: 
#          p_y    : approximate conditional density
p_y = function(theta, prior, K, data) {
        
    y = data$y    # (1 x 1) response variable
    x = data$x    # (D x 1) covariate vector 

    p_y = 0
    
    for (k in 1:K) {
        # k-th gaussion: N ( y | x'beta_k, tau_k^{-1} )
        tau_k_inv  = theta$b_k[k] / theta$a_k[k]    # precision component
        beta_k     = theta$m_k[,k]                  # coefficient vector
        mu_k       = c(t(x) %*% beta_k)             # mean component
        p_y        = p_y + theta$pi_k[k] * dnorm(y, mean = mu_k, tau_k_inv)
    }

    return(p_y)
    
} # end approxDensity() function


# approx_p_y() : evaluate each of the y values using the mixture density 
# input: 
#          y_seq    :
#          x        : 
#          theta    : 
#          prior    : 
# output: 
#          p_y_vec  :  p_y(y) values for each of the y_seq inputs 
approx_p_y = function(y_seq, x, theta, prior) {





} # end approx_p_y() function




