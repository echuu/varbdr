
# varbdr.R -- bayesian density regression with covariate independent weights,
# sparsity assumption the gaussian components of the mixture of normals

varbdr = function(prior, theta) {
    
    
    for (i in 2:theta$max_iter) {
        
        # variational update
        theta = updateQ(prior, theta)            
        
        if (converge(prior, theta)) {
            break
        }    
        
        theta$curr = theta$curr + 1
        
    } # end of CAVI loop
    
    
    return(theta)
    
} # end of varbdr() function


