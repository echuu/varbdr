
## beta2.R -- generate data of the form:

##     X_n = (X1, ... , X4) ~ Unif[0.05, 0.95]^4
##     y_n ~ Beta(5 * X2 * exp(2 * X1), 5 * X3 + 3 * X4)

# beta1() 
# input:  
#          N     : number of observations/data points to generate
#          seed  : seed for RNG 
# output: 
#          X     : (N x 2) design matrix, covariates stored row-wise
#          y     : (N x 1) vector of responses
#          shape : (N x 2) shape parameters for each y_n

beta2 = function(N, seed = 100) {
    
    set.seed(seed)
    D = 4                   # dimension of covariates
    y = numeric(N)          # (N x 1) vector of response variables
    X = matrix(0, N, D)     # (N x D) covariates for y_n are stored row-wise  
    
    
    shape = matrix(0, N, 2) # shape parameters for beta distribution
    
    for (n in 1:N) {
        # X_n = (X_n1, X_n2, ... , X_nD) ~ Unif [0.05, 0.95]^D
        X[n,] = runif(D, 0.05, 0.95)
    }
    
    # store shape parameters for beta distribution
    shape[,1] = 5 * X[,2] * exp(2 * X[,1])
    shape[,2] = 5 * X[,3] + 3 * X[,4]
        
    # y_n | x_n ~ Beta(shape1, shape2)
    y = rbeta(N, shape1 = shape[,1], shape2 = shape[,2])
     
    
    
    return(list(X = X, y = y, shape = shape))
    
}