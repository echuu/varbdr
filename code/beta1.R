
## beta1.R -- generate data of the form:

##     X_n = (X1 , X2) ~ Unif[0.05, 0.95]^2
##     y_n ~ Beta(4 * X1 + 3 * X2^2, 10 * X2)

# beta1() 
# input:  
#          N     : number of observations/data points to generate
#          seed  : seed for RNG 
# output: 
#          X     : (N x 2) design matrix, covariates stored row-wise
#          y     : (N x 1) vector of responses
#          shape : (N x 2) shape parameters for each y_n

beta1 = function(N, seed = 100) {
    
    set.seed(seed)
    D = 2                   # dimension of covariates
    y = numeric(N)          # (N x 1) vector of response variables
    X = matrix(0, N, D)     # (N x D) covariates for y_n are stored row-wise  
    
    shape = matrix(0, N, 2) # shape parameters for beta distribution 
    
    for (n in 1:N) {
        # X_n = (X_n1, X_n2, ... , X_nD) ~ Unif [0.05, 0.95]^D
        X[n,] = runif(D, 0.05, 0.95) 
    }
     
    # store shape parameters for beta distribution
    shape[,1] = 4 * X[,1] + 3 * X[,2]^2
    shape[,2] = 10 * X[,2]
        
    # y_n | x_n ~ Beta(shape1, shape2)
    y = rbeta(N, shape1 = shape[,1], shape2 = shape[,2])
     
    return(list(X = X, y = y, shape = shape))
    
} # end beta1() function


old_beta1 = function(N, seed = 100) {
    
    set.seed(seed)
    D = 2                   # dimension of covariates
    y = numeric(N)          # (N x 1) vector of response variables
    X = matrix(0, N, D)     # (N x D) covariates for y_n are stored row-wise  
    
    shape = matrix(0, N, 2) # shape parameters for beta distribution 
    
    for (n in 1:N) {
        # X_n = (X_n1, X_n2, ... , X_nD) ~ Unif [0.05, 0.95]^D
        X[n,] = runif(D, 0.05, 0.95) 
        
        # store shape parameters for beta distribution
        shape[n, 1] = 4 * X[n, 1] + 3 * X[n, 2]^2
        shape[n, 2] = 10 * X[n, 2]
        
        # y_n | x_n ~ Beta(shape1, shape2)
        y[n] = rbeta(1, shape1 = shape[n,1], shape2 = shape[n, 2])
    } 
    
    
    return(list(X = X, y = y, shape = shape))
    
} # end beta1() function




# end beta1.R
