
## approxDensity.R -- for ALL cases
## generate the approximating density using a mixture of experts, using
## the model parameters from CAVI

## functions included:
##   (1) py_0          : approximate mixture density (cov-indep)
##   (2) py_bouch      : approximate mixture density (cov-dep, bouchard)
##   (3) plotDensities : overlay true density and approximate mixture density


library(ggplot2)
library(reshape2)

# py_0(): mixture density used to approximate the true density,
#        covariate-INDEPENDENT weights
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
py_0 = function(theta, K, data) {
    
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
    
} # end py_0() function



# py_bouch(): mixture density used to approximate the true density,
#             using covariate-dependent weights (Bouchard bound)
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
py_bouch = function(theta, K, data) {
    
    y = data$y    # (1 x 1) response variable
    x = data$x    # (D x 1) covariate vector 
    
    
    # calculate pi_k
    x_gamma = t(x) %*% theta$gamma_k           # (1 x K)
    pi_k = exp(x_gamma - log_sum_exp(x_gamma)) # (1 x K) : normalize the weights 
    
    p_y = 0
    for (k in 1:K) {
        # k-th gaussion: N ( y | x'beta_k, tau_k^{-1} )
        tau_k_inv  = 1 / theta$tau_k[k]             # precision component
        beta_k     = theta$beta[,k]                 # coefficient vector
        mu_k       = c(t(x) %*% beta_k)             # mean component
        p_y        = p_y + pi_k[k] * dnorm(y, mean = mu_k, sd = sqrt(tau_k_inv))
    }
    
    return(p_y)
    
} # end py_bouch() function



# compareDensities() : compare approximating densities to true (beta1) density
# input:  
#          y_grid     :  sequence of y-values evaluated using each density
#          x          :  covariates for the y (response) values
#          true_d     :  true density function
#          params     :  parameters for the true density function
#          den_list   :  list of each approximating density to use to evaluate
#          true_den   :  true density
#          den_label  :  names of each density (plotting purposes)
#          theta_list :  variational parameters
#          K          :  number of clusters used in the mixture density
#        
# output: 
#          d_plot     :  ggplot of overlayed densities
compareDensities = function(y_grid, x, 
                            true_d, params,
                            den_list, den_label, 
                            theta_list, K) {
    
    N_evals    = length(y_grid)           # number of evaluations
    n_den      = length(den_label)        # number of densities to evaluate
    data_ygrid = list(y = y_grid, x = x)  # fixed for all evaluations
    
    # approx_mat: store density evaluation results (N_evals) x (n_den + 1)
    #     store y_grid in the first column
    #     store the density evaluations in rows 2 to n_den + 1
    #     extra column for the true density evaluations
    approx_mat = matrix(0, nrow = N_evals, ncol = n_den + 1)
    # approx_mat[,n_den] = y_grid
    
    # evaluate using each of the approximating densities
    for (d_i in 1:n_den) {
        approx_mat[, d_i] = den_list[[d_i]](theta_list[[d_i]], K, data_ygrid)
    }
    
    # true density -- this line needs to be changed depending on 
    # what the true density is
    approx_mat[,n_den + 1] = true_d(y_grid, 
                                    shape1 = params[[1]], shape2 = params[[2]])
    
    approx_df = data.frame(y_grid, approx_mat)
    names(approx_df) = c("y", den_label, "true")
    
    approx_df_long = melt(approx_df, measure.vars = c(den_label, "true"))
    
    ggplot(approx_df_long, aes(x = y, y = value, colour = variable)) + 
        geom_point(size = 0.7) + labs(x = "y", y = "p(y)") + theme_bw()
    
}

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
                         approx_d, theta, K) { 
    
    data_ygrid = list(y = y_grid, x = x)  # list req'd to use approx density fcn
    N_evals    = length(y_grid)           # number of evaluations
    
    # evaluate points using approx density
    p_ygrid    = approx_d(theta, K, data_ygrid)
    approx_df  = data.frame(x = y_grid, y = p_ygrid)   # store y, f(y) values
    
    # true density plot (indexed by n in the shape parameters)
    beta_n = stat_function(aes(x = y_grid, y = ..y..), 
                           fun = dbeta, colour = 'black', n = N_evals, size = 1,
                           args = params)  
    
    
    p_line = geom_line(aes(x = approx_df$x, y = approx_df$y),
                       colour = "blue", size = 0.9)
    # p = ggplot(approx_df, aes(x, y)) + geom_point(colour = 'blue', size = 0.9)
    p = ggplot(approx_df, aes(x, y)) + 
        geom_line(aes(x = approx_df$x, y = approx_df$y),
                  colour = "blue", size = 0.9)
    p = p + beta_n
    
    return(list(overlay = p, approx = p_line))
} # end of plotDensities() function



# densityCurve(): return density curve for each set of N covariates
# input:  
#          approx_d :  approx density function
#          theta    :  variational parameters
#          X        :  (N x D) design mat w/ covariates stacked row-wise
#          K        :  number of clusters used in the mixture density
#          y_grid   :  sequence of y-values evaluated using true/approx density
# output: 
#          approx   :  N x len dataframe of p_y evaluations -> density curve
densityCurve = function(approx_d, theta, X, K,
                        y_grid = seq(0, 1, len = 500)) {
    
    N = nrow(X)
    
    approx_t = matrix(0, nrow = N, ncol = length(y_grid))
    
    # obtain density curve estimates for each of the n observations in 
    # design matrix X
    for (n in 1:N) {
        data_ygrid  = list(y = y_grid, x = X[n,])
        # print(data_ygrid$x)
        p_ygrid     = approx_d(theta, K, data_ygrid)
        # approx_df   = data.frame(x = y_grid, y = p_ygrid)
        approx_t[n,] = p_ygrid
    }
    
    return(approx_t)
} # end of densityCurve() function



# queryDensity() : for a given observation(s), generate the corresponding 
#                  density plot(s) and overlay with the true density
# input:  
#          X        :  (N x D) design mat w/ covariates stacked row-wise
#          n_set    :  row indices of X to plot conditional densities
#          params   :  matrix of parameters for the true density
#          true_d   :  true density
#          approx_d :  approx density function
#          d_label  :  vector of strings of names of each of the algorithms
#          theta    :  list of variational parameters (for diff vb algos)
#          K        :  number of clusters used in the mixture density
#          y_grid   :  sequence of y-values evaluated using true/approx density
# output: 
#          approx   :  N x len dataframe of p_y evaluations -> density curve
queryDensity = function(X, n_set, params, true_d,
                        approx_d, d_label, theta, K,
                        y_grid = seq(0, 1, len = 500)) {
    
    num_plots = length(n_set)
    plot_list = vector("list", num_plots)
    
    for (i in 1:num_plots) {
        
        param_n = list(shape1 = params[n_set[i], 1], 
                       shape2 = params[n_set[i], 2])
        
        plot_list[[i]] = compareDensities(y_grid, X[n_set[i],], 
                                          true_d, param_n,
                                          approx_d, d_label, 
                                          theta, K)
    }
    
    return(plot_list)
}



# multiplot() : plot multiple ggplot2 objects in the same figure
multiplot = function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots = c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots == 1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get i,j matrix positions of the regions that contain this subplot
            matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
} # end multiplot() function





# end of density.R
