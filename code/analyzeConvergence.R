




plotParams = function(theta, K, D) {
    
    iter = theta$curr - 1 # num iters to reach convergence
    
    V_k = t(theta$V_k_hist[,1:iter]) # iter x K
    Q_k = t(theta$Q_k_hist[,1:iter]) # iter x K
    
    V_k_long = melt(V_k)
    Q_k_long = melt(Q_k)
    
    names(V_k_long) = c("iter", "cluster", "value")
    names(Q_k_long) = c("iter", "cluster", "value")
    
    V_k_plot = ggplot(V_k_long, aes(x = iter, y = value, 
                                    col = factor(cluster))) + geom_point()  + 
        ggtitle("Frob norm for V_k vs. Iter") + theme_bw()
    
    Q_k_plot = ggplot(Q_k_long, aes(x = iter, y = value, 
                                    col = factor(cluster))) + geom_point()  + 
        ggtitle("Frob norm for Q_k vs. Iter") + theme_bw()
    
    
    m_k_plot = vector('list', K)
    mu_k_plot = vector('list', K)
    for (k in 1:K) {
        
        m_k  = theta$m_k_hist[,,k][1:iter]   # (D x iter) matrix
        mu_k = theta$mu_k_hist[,,k][1:iter]  # (D x iter) matrix
        
        if (D > 1) {
            
            m_k_long = melt(m_k)    # each column will be its own id
            mu_k_long = melt(mu_k)
            
            m_k_plot[[k]] = ggplot(m_k_long, aes(x = 1:iter, y = value,
                                                 col = factor(variable))) +
                geom_point() + ggtitle(paste("Value of m_", k, sep = '')) + 
                theme_bw()
            
            mu_k_plot[[k]] = ggplot(mu_k_long, aes(x = 1:iter, y = value,
                                                   col = factor(variable))) +
                geom_point() + ggtitle(paste("Value of mu_", k, sep = '')) + 
                theme_bw()
            
        } else {
            m_k_df = data.frame(i = 1:iter, m_k = m_k)   
            mu_k_df = data.frame(i = 1:iter, mu_k = mu_k)  
            
            m_k_plot[[k]] = ggplot(m_k_df, aes(i, m_k)) + geom_point() + 
                ggtitle(paste("Value of m_", k, sep = '')) + theme_bw()
            
            mu_k_plot[[k]] = ggplot(mu_k_df, aes(i, mu_k)) + geom_point() + 
                ggtitle(paste("Value of mu_", k, sep = '')) + theme_bw()
        }
        
        
    }
    
    
    plot_list = list(m_k_plot = m_k_plot, mu_k_plot = mu_k_plot, 
                     V_k_plot = V_k_plot, Q_k_plot = Q_k_plot)
    return(plot_list)
}





