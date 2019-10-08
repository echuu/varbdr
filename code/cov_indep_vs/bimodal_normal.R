

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
d = data$d

xy_df = data.frame(x = X[,d], y = y)
ggplot(xy_df, aes(x, y)) + geom_point()

# plot conditional densities for different values of x
y_grid   = seq(-3, 2, length.out = 5000)
num_pred = 3 # num of conditional densities to be evaulated
x        = matrix(runif(num_pred * D), num_pred, D)
x[,d]    = c(0, 0.5, 0.9) # fill in d-th column with the non-noise covariate
f_yx     = matrix(0, nrow = length(y_grid), ncol = num_pred) # (1000 x 3)


yx_curves = matrix(0, nrow = length(y_grid), ncol = num_pred) # (1000 x 3)
yx_plots  = list()

for(i in 1:num_pred) {
    # compute the conditional density for each (y, x) pair using the 
    # mixture of normals density 
    # write this in the simulation.R file and generalize the probability p 
    # used in the gmm_sim_xy() function
    
    # for each value of x, find the conditional density and plot its curve
    # y_eval = data.frame(y = y_grid, x = x[i])
    yx_curves[,i] = gmm_pdf(y_grid, x[i,]) # gmm_pdf() in simulation.R file
    
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
source('misc.R')   # load convergence functions


## initialize model hyper-parameters
alpha_0 = rep(1 / K, K);                         # pi_k params
m_0 = 0; xi_0 = 10;                              # beta params             
pi_d = 0.3;                                      # prior prob of inclusion
a_0 = 0.5; b_0 = 0.5;                                # tau params
VERBOSE = TRUE; tol = 1e-3;  max_iter = 50;      # algo/converge params


## initialize variational parameter objects
prior = initPriors(y, X, K, 
                   alpha_0,                    # pi   : concentration param
                   m_0, xi_0,                  # beta : mean, scale for slab
                   pi_d,                       # beta : prior pip
                   a_0, b_0)                   # tau  : shape, rate

# theta = initVarParams_indep(y, X, N, D, K, max_iter = max_iter)

source('updateFunctions.R') # update functions for variational distributions
source('updateQ.R')         # performs one update for each distribution
source('varbdr.R')          # main cavi loop
source('elbo.R')

# run the cavi algorithm
theta = initVarParams_indep(y, X, N, D, K, max_iter = max_iter)
theta_star = varbdr(prior, theta)

# theta_star stores variational parameters after convergence

## TODO: -----------------------------------------------------------------------

# obtain the approximate conditional density using the x-values passed into
# mixture of experts model

source('approxDensity.R')

# store the conditional density approximations (col-wise per x input)
py_x = matrix(0, nrow = length(y_grid), ncol = num_pred) # (1000 x 3)

for (i in 1:num_pred) {
    # evaluate each of the conditional densities
    py_x[,i] = approx_f_xy(prior, theta_star, y_grid, x[i,])
}

# plot the approx conditional densities overlayed with the 'true' condititional
# densities 

cd_plots = vector("list", num_pred)

for (i in 1:num_pred) {
    # data frame for y, p(y|x) evaluations, true cond. densities f(y|x)
    df_y_py = data.frame(y = y_grid, py = py_x[,i], fy = yx_curves[,i])
    
    # create conditional density plot
    cd_plots[[i]] = ggplot(df_y_py, aes(x = y, y = py)) + 
        geom_point(size = 0.5)
}


# overlay the true density (blue) on the approximate density (black)
cd_plots[[1]] + geom_line(aes(x = y_grid, y = yx_curves[,1]), color = 'blue')

cd_plots[[2]] + geom_line(aes(x = y_grid, y = yx_curves[,2]), color = 'blue')

cd_plots[[3]] + geom_line(aes(x = y_grid, y = yx_curves[,3]), color = 'blue')

multiplot(plotlist = cd_plots, cols = num_pred)


# numerical figure to quantify discrepancy

# uncertainty estimate (?) need to look at how C&S + others do this

# quantify accuracy of the variable selection


## DONE: -----------------------------------------------------------------------







