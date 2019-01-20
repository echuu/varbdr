
## misc.R

## misc calculations performed in CAVI for covariate-INDEPENDENT case

# computes the log of the normalizing constant of the Dirichlet density
# log(1/B(alpha)) = log(gamma(sum(alpha_k)) / prod(gamma(alpha_k)))
# intput:
#         alpha_k : (K x 1) or (1 x 1) if symm. Dirichlet; concentration paramt
#         K       : (1 x 1) dimension of dirichlet random variable
# output: 
#         log of the normalizing constant of Dirichlet density
log_dir_const = function(alpha_k, K) {
  
    log_const = 0
    
    if (K != length(alpha_k)) {      # symmetric dirichlet, alpha_k is scalar
        alpha_k = rep(alpha_k, K)    # in this case, turn it into k-dim vec
    }
    
    log_const = lgamma(sum(alpha_k)) - sum(lgamma(alpha_k))
      
    return(log_const)
    
} # end log_dir_const() function




