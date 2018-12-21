
# testResults.R
run = function() {
    set.seed(12)
    
    # setwd to source file location
    setwd("C:/Users/chuu/varbdr/cov_independent_code/r_implementation")
    source("displayResults.R")
    source("gmVB_0.R")
    
    
    # View(faithful) # 272 x 2 : duration of eruption, waiting time b/w eruptions
    
    X = as.matrix(faithful)
    K = 25        # Number of clusters
    
    K = 25
    alpha_0 = 1e-5
    m_0 = c(colMeans(X))
    beta_0 = 1 
    nu_0 = NCOL(X) + 50
    W_0 = diag(100, NCOL(X))
    max_iter = 1001
    epsilon_conv = 1e-4
    is_animation = TRUE
    VERBOSE = TRUE
    
    
    
    
    # Run vb-gmm model model
    vb_gmm_model = vb_gmm(X = X, K = K, alpha_0 = 1e-5, max_iter = 1001, 
                           is_animation = TRUE, VERBOSE = TRUE)
}




data.grid = expand.grid(x = seq(from = min(X[,1]) - 2, 
                                 to = max(X[,1]) + 2, length.out = 100), 
                         y = seq(from = min(X[,2]) - 8, 
                                 to = max(X[,2]) + 2, length.out = 100))

# generate from predictive desntiy
q.samp = cbind(data.grid, z = mixture_pdf_t(vb_gmm_model,data.grid))

ggplot() + 
    geom_point(data = data.frame(X), mapping = aes(eruptions, waiting)) + 
    geom_contour(data = q.samp, 
                 mapping = aes(x = x,y = y, z = z,
                               colour = ..level..), binwidth = 0.001) + 
    gg_theme()
