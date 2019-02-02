
## approxDensity.R -- covariate INDEPENDENT case
## generate the approximating density using a mixture of experts, using
## the model parameters from CAVI

## functions included:
##   (1) p_y           : approximate mixture density
##   (2) plotDensities : overlay true density and approximate mixture density


library(ggplot2)

# p_y(): mixture density used to approximate the true density
# input:  
#          theta  : variational parameters
#          prior  : prior-related variables
#          K      : number of clusters used in the mixture density
#          data   : LIST containing the response (y) and the covariates (x)
#                   note: y is a (N x 1) vector, x is (D x 1) vector
#                         every y value will be evaluated using the same
#                         covariate vector x
# output: 
#          p_y    : approximate conditional density
p_y = function(theta, prior, K, data) {
        
    y = data$y    # (1 x 1) response variable
    x = data$x    # (D x 1) covariate vector 

    p_y = 0
    
    for (k in 1:K) {
        # k-th gaussion: N ( y | x'beta_k, tau_k^{-1} )
        tau_k_inv  = 1 / theta$tau_k[k]             # precision component
        beta_k     = theta$beta[,k]                 # coefficient vector
        mu_k       = c(t(x) %*% beta_k)             # mean component
        p_y        = p_y + theta$pi_k[k] * dnorm(y, mu_k, sqrt(tau_k_inv))
    }

    return(p_y)
    
} # end approxDensity() function



# plotDensities(): overlay the true density and the (approx) mixture density
# input:  
#          y_grid   :  sequence of y-values evaluated using true/approx density
#          x        :  covariates for the y (response) values
#          true_d   :  true density function
#          params   :  parameters for the true density function
#          approx_d :  approx density function
#          theta    :  variational parameters
#          prior    :  prior-related variables
#          K        :  number of clusters used in the mixture density
# output: 
#          d_plot   :  ggplot of overlayed densities (true: red, approx: blue)
plotDensities = function(y_grid, x,
                         true_d, params,
                         approx_d, theta, prior, K) { 
    
    data_ygrid = list(y = y_grid, x = x)  # list req'd to use approx density fcn
    N_evals    = length(y_grid)           # number of evaluations
    
    # evaluate points using approx density
    p_ygrid    = approx_d(theta, prior, K, data_ygrid)
    approx_df  = data.frame(x = y_grid, y = p_ygrid)   # store y, f(y) values
    
    # true density plot (indexed by n in the shape parameters)
    beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                           fun = dbeta, colour = 'black', n = N_evals, size = 1,
                           args = params)  
    
    
    p = ggplot(approx_df, aes(x, y)) + geom_point(colour = 'blue', size = 0.9)
    p = p + beta_n
    
    return(p)
}




densityCurve = function(approx_d, theta, prior, x, K,
                        y_grid = seq(0, 1, len = 500)) {
    
    data_ygrid = list(y = y_grid, x = x)  
    p_ygrid    = approx_d(theta, prior, K, data_ygrid)
    approx_df  = data.frame(x = y_grid, y = p_ygrid)
    
    approx_i   = geom_line(aes(x = approx_df$x, y = approx_df$y),
                           colour = "blue", size = 0.9)
    
    return(approx_i)
    
}









# end of approxDensity.R file
