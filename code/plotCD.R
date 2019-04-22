



generatePlots = function(py_df, true_den, k_den) {
    
    plot_list = vector('list', length(py_df))
    
    py_names = c("vb_d")
    if (!is.null(true_den)) {
        py_names = c(py_names, "true_d")
    }
    if (!is.null(k_den)) {
        py_names = c(py_names, "kern_d")
    }
    
    
    for (i in 1:length(py_df)) {
        
        df_i_long = melt(py_df[[i]], 
                         measure.vars = py_names)
        
        
        if (!is.null(true_den)) {
            df_i_long = df_i_long %>% mutate(is_approx = 'true')
            df_i_long$is_approx[df_i_long$variable != 'true_d'] = 'approx'
        }
        
        
        plot_list[[i]] = ggplot(df_i_long, 
                                aes(x = y, y = value, col = variable)) + 
            geom_line(size = 0.8, 
                      aes(linetype = factor(df_i_long$is_approx))) +
            labs(x = "y", y = "p(y)") + theme_bw() +
            theme(legend.position = "none")
        
    }
    
    
    return(plot_list)
    
    
}


plotCD = function(theta, K, x_in, y_grid, true_den = NULL, k_den = NULL) {
    
    
    # if intercept is fitted, then need to include a 1 before each of the 
    # input x values
    if (theta$intercept) {    
        x_mat = as.matrix(cbind(1, x_in)) 
    } else {    
        x_mat = as.matrix(x_in)              
    } 
    
    N_evals    = length(y_grid)             # number of function evals for f(y)
    
    num_plots = length(x_in)
    py_df     = vector("list", num_plots)   # store each of the dataframes
    plots_py  = vector("list", num_plots)   # store each of the resulting plots
    
    # generate density plots for each x input : for each x, we compute the
    # density y along y_grid; this gives us f(y | x)
    
    for (i in 1:length(x_in)) {
        
        py_vb   = NULL
        py_true = NULL
        py_kern = NULL
        
        x0 = x_mat[i,] # pass in i-th covariate - depending on intercept, this
                       # will look like either: (1) [x] or (2) [1, x]
        
        # compute density evalulations f ( y | x[i] )
        
        # (0) compute vb approximation
        xy_list = list(y = y_grid, x = x0)   # fixed x for all evaluations
        py_vb = py_bouch(theta, K, xy_list)
        
        
        if (length(x0) > 1) {     # omit the intercept term when evaluating TRUE
            x0 = x0[2:length(x0)] # density bc true density only requires x
        }
        
        # (1) compute true density (if applicable)
        if (!is.null(true_den)) {
            py_true = true_den(y_grid, x0)
        }
        
        # (2) compute kernel density (if applicable)
        if (!is.null(k_den)) {
            y_eval  = data.frame(y = y_grid, x = x0) # predict requires df
            py_kern =  predict(k_den, newdata = y_eval)
        } 
        
        
        # store the values for p(y) (between 1 and 3 columns) in a matrix
        py_df[[i]] = data.frame(y = y_grid, vb_d = py_vb)
        
        py_df[[i]]$true_d = py_true
        py_df[[i]]$kern_d = py_kern
        # true_d = py_true, kern_d = py_kern
        
        # approx_df_long = melt(py_df[[i]], measure.vars = c(den_label, "true"))
        
    }
    
    
    cd_plots = generatePlots(py_df, true_den, k_den)
    
    
    return(list(py_df = py_df, cd_plots = cd_plots))
    
} # end of plotCD() function





