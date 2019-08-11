
# varbdr.R -- bayesian density regression with covariate independent weights,
# sparsity assumption the gaussian components of the mixture of normals

varbdr = function(prior, theta) {
    
    conv_status = ''
    
    for (i in 2:theta$max_iter) {
        
        # variational update
        theta = updateQ(prior, theta)            
        
        # compute variational lower bound
        theta = computeELBO(prior, theta)
        
        if (theta$converge) {
            conv_status = 'elbo converged'
            break
        }    
        
        # convergence status messages -- set VERBOSE = FALSE to suppress
        if (theta$VERBOSE) {
            print(paste('iter:', theta$curr - 1, conv_status, sep = ' '))
        }
        
        theta$curr = theta$curr + 1
        
    } # end of CAVI loop
    
    
    return(theta)
    
} # end of varbdr() function


