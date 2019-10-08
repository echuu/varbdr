
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
            
            # include the elbo plot -- move this into elbo file later
            
            elbo_df = data.frame(iter = 1:(i-1), elbo = theta$L[2:i])
            
            theta$elbo_plot = ggplot(elbo_df, aes(iter, elbo)) + geom_point() +
                scale_x_discrete(limits = c(1:(i-1))) + 
                theme_bw()
            
            print(paste('iter:', theta$curr - 1, conv_status, 
                        '--------------------------------', sep = ' '))
            break
        }    
        
        # convergence status messages -- set VERBOSE = FALSE to suppress
        if (theta$VERBOSE) {
            print(paste('iter:', theta$curr - 1, conv_status, 
                        '--------------------------------', sep = ' '))
        }
        
        theta$curr = theta$curr + 1
        
    } # end of CAVI loop
    
    
    return(theta)
    
} # end of varbdr() function


