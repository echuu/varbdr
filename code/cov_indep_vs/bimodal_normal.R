

source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
# source(DP_BDR)    # so that we have access to data

setwd("~/varbdr/code/cov_indep_vs")
source("simulation.R") # generate data from conditional density
source("~/varbdr/code/density.R")

# generate data
N = 1000   # number of observations
D = 6      # dimension of covariates
K = 2      # number of clusters
data = gmm_sim_xy(N, D)
y = data$y
X = data$X

xy_df = data.frame(x = X[,1], y = y)
ggplot(xy_df, aes(x, y)) + geom_point()

# plot conditional densities for different values of x
y_grid   = seq(-3, 2, length.out = 1000)
x        = c(0, 0.5, 0.9)
f_yx     = matrix(0, nrow = length(y_grid), ncol = length(x)) # (1000 x 3)

yx_curves = matrix(0, nrow = length(y_grid), ncol = length(x))
yx_plots = list()

for(i in 1:length(x)) {
    # compute the conditional density for each (y, x) pair using the 
    # mixture of normals density 
    # write this in the simulation.R file and generalize the probability p 
    # used in the gmm_sim_xy() function
    
    # for each value of x, find the conditional density and plot its curve
    # y_eval = data.frame(y = y_grid, x = x[i])
    yx_curves[,i] = gmm_pdf(y_grid, x[i]) # gmm_pdf() in simulation.R file
    
    # create ggplot() object and store in yx_plots'
    df_fy = data.frame(y = y_grid, f_yx = yx_curves[,i])
    yx_plots[[i]] = ggplot(df_fy, aes(y, f_yx)) + geom_line(size = 0.8) + 
        labs(x = 'y', y = '') + theme_bw()
    
}

# plot the 3 conditional density curves
multiplot(plotlist = yx_plots, cols = 3)

# ------------------------------------------------------------------------------

source("initPr.R") # load initPriors() function
source("initVP.R") # load initVarParams_indep() function

## initialize model hyper-parameters
alpha_0 = rep(1 / K, K);                         # pi_k params
m_0 = 0; xi_0 = 0.3;                             # beta params             
pi_d = 0.3;                                      # prior prob of inclusion
a_0 = 1; b_0 = 1;                                # tau params
VERBOSE = TRUE; tol = 1e-3;  max_iter = 50;      # algo/converge params


## initialize variational parameter objects

prior = initPriors(y, X, K, 
                   alpha_0,                    # pi   : concentration param
                   m_0, xi_0,                  # beta : mean, scale for slab
                   pi_d,                       # beta : prior pip
                   a_0, b_0)                   # tau  : shape, rate

theta = initVarParams_indep(y, X, N, D, K)

source('updateFunctions.R') # update functions for variational distributions
source('updateQ.R')         # performs one update for each distribution
source('misc.R')            # load convergence functions
source('varbdr.R')          # main cavi loop

# run the cavi algorithm
theta = initVarParams_indep(y, X, N, D, K)
theta_star = varbdr(prior, theta)


