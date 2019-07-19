
# varbdr_vs.R -- bayesian density regression with covariate dependent weights,
# sparsity assumption the gaussian components of the mixture of normals


# in the prev implementation, prior and theta were both created inside the
# varbdr function, but for now, we initialize those outside so that debugging 
# is easier (don't have to re-initialize prior, theta objects to perform the 
# actual inference procedure)
varbdr = function(prior, theta) {
    
    for (i in 2:theta$max_iter) {
        
        # TODO: implement eStep_vs() function
        theta = eStep_vs(prior, theta)
        theta = mStep_vs(prior, theta)
        
        theta$L = elbo_vs(prior, theta)
        
        # TODO: implement converge() -- previously called checkELBO()
        if (converge(prior, theta)) {
            break
        }
        
        theta$curr = theta$curr + 1
        
    } # end of CAVI loop
    
    
    # return the set of updated parameters
    return(theta)
    
} # end of varbdr() function







