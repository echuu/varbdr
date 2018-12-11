

## displayResults.R


## functions that display results

# Define ggplot2 theme
gg_theme <- function(){
    p <- theme(
        plot.title = element_text(size = 20,face = 'bold',
                                  margin = margin(0,0,3,0), hjust = 0.5),
        axis.text = element_text(size = rel(1.05), color = 'black'),
        axis.title = element_text(size = rel(1.45), color = 'black'),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
        axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1.4, 'lines'),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.major = element_line(colour = "gainsboro"),
        panel.background = element_blank()
    )
    return(p)
}

# end of displayResults.R


## various checks on the ELBO
## input:
    # VERBOSE      : logical, if TRUE, then progress printed each iter
    # i            : current iteration
    # max_iter     : max # of iterations
    # epsilon_conv : tolerance used to assess convergence
checkELBO = function(VERBOSE, i, max_iter, L, epsilon_conv) {
    
    CONVERGE = FALSE
    
    # display iteration, ELBO, change in ELBO
    if (VERBOSE) {
        cat("It:\t",i,"\tLB:\t",L[i], "\tLB_diff:\t",L[i] - L[i - 1],"\n")
    }
    
    # Check if lower bound decreases
    if (L[i] < L[i - 1]) { 
        message("Warning: Lower bound decreases!\n")
    }
    # Check for convergence
    if (abs(L[i] - L[i - 1]) < epsilon_conv) { 
        CONVERGE = TRUE 
    }
    
    # Check if VB converged in the given maximum iterations
    if (i == max_iter) {
        warning("VB did not converge!\n")
    }
    
    return(CONVERGE)
} # end of checkELBO() function







