

## eStep_vs.R -- covariate dependent with sparsity assumption
## perform one iteration of the variational e-step



# input: 
#          theta : list of variational parameters
# output:
#          theta : list of variational parameters with 
#                  variational parameters updated

eStep_vs = function(prior, theta) {
    
    # (1) update: q(Z) = prod_n prod_k q(z_nk)
    theta = rnkUpdate(prior, theta)
    
    return(theta)
    
} # end of eStep_vs() function
