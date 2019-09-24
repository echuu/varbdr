
# testFunctions.R -- test functions for SSCI-MoE


source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR)    # so that we have access to data

setwd("~/varbdr/code/cov_indep_vs")
source("misc.R")



# generate data
N = 500
D = 1
K = 2
synth_data_1d = r_dpmix2(N)
y = synth_data_1d$y
X = synth_data_1d$X
# X = scale(X, scale = TRUE)

intercept = FALSE;     # dated setting? used to check that 1st column is 1 vec

# initialize parameters needed to pass into initPriors() function
# TODO: figure out the appropriate value to initialize pi_d with
# C&S implementation seems to normalize the probabilities over all features (?)

alpha_0 = rep(1 / K, K);                         # pi_k params
m_0 = 0; xi_0 = 0.3;                             # beta params             
pi_d = 0.3;                                      # prior prob of inclusion
a_0 = 1; b_0 = 1;                                # tau params
VERBOSE = TRUE; tol = 1e-3;  max_iter = 1e4;     # algo/converge params


# test initPriors() function
source("initPr.R")
prior = initPriors(y, X, K, 
                   alpha_0,                    # pi   : concentration param
                   m_0, xi_0,                  # beta : mean, scale for slab
                   pi_d,                       # beta : prior pip
                   a_0, b_0)                   # tau  : shape, rate
                   

# test initVarParams() function


# check these unconditional exp/cov quantities
theta$Sigma_k
theta$m_k
theta$var_beta_d

source('updateFunctions.R')

# (0) test rnkUpdate() function
theta = rnkUpdate(prior, theta)

# (1) test ssUpdate() function
theta = ssUpdate(prior, theta)



## ---- testing ------------------------------------------------------------- ##

source("initVP.R")
source('updateFunctions.R')
theta = initVarParams_indep(y, X, N, D, K)
theta = rnkUpdate(prior, theta)
theta = ssUpdate(prior, theta)

# (1.1) test betak_update() function
theta = betak_update(prior, theta)

# check these unconditional exp/cov quantities
theta$Sigma_k
theta$m_k
theta$var_beta_d

## 8/11 progress update --- 

## TODO: betak_update() needs to be called before the first e-step update
## problem -- betak_upate() requires prior, theta object, but the update would
## be called from within initVP() -- solution (?) : call betak update the rest
## of the updates since it's only used "after" the spike and slab update which
## is the last step of the coordinate ascent updates

## solution: put betak_update() code into the initVP() so that all the values 
## are updated accordingly. Then call betak_update() as planned after after
## the ssUpdate()

## other issues: lambda_d = 1 on the first iteration causing errors in the elbo
## computatio b/c there is log(lambda_d - 1) = log(0) = NaN


## ---- testing ------------------------------------------------------------- ##


# (2) test precUpdate() function
theta = precUpdate(prior, theta)

# (3) test wtUpdate() function
theta = wtUpdate(prior, theta)


# ---------------------------------------------------------------------------- #

# test updateQ() function
source('updateQ.R')
theta_p = updateQ(prior, theta)

# TODO: update m_k after updating m_d; this is not just the transpose of m_d
# since m_d = E(beta_d | omega_d = 1), whereas we're interpreting m_k as
# E(beta_k), the unconditional expectation (which is used in updates for
# q(tau), q(z)), nor used in the ELBO, E[ln p(y | -)] written in terms fo beta_d

# TODO: update Sigma_k (where did I do this before?)

### takeaway: m_k, Sigma_k need to be updated in the same chunk of the code 
###           as the other spike and slab parameters since they are used 
###           in successive calculations (next iteration)

source('elbo.R')
theta_elbo = computeELBO(prior, theta)


## test varbdr() function
source('varbdr.R')





# ---------------------------------------------------------------------------- #


source("~/varbdr/code/globals.R")
setwd(HOME_DIR)
source(DP_BDR)    # so that we have access to data

setwd("~/varbdr/code/cov_indep_vs")
source("misc.R")



# generate data from conditional density
source("simulation.R")

# generate data
N = 1000   # number of observations
D = 4      # dimension of covariates
K = 3      # number of clusters
data = gmm_sim_xy(N, D)
y = data$y
X = data$X
# X = scale(X, scale = TRUE)

# intercept = FALSE;     # dated setting? used to check that 1st column is 1 vec

# plot y ~ x
xy_df = data.frame(x = X[,1], y = y)
ggplot(xy_df, aes(x, y)) + geom_point()

# plot conditional densities for different values of x
y_grid   = seq(-3, 2, length.out = 1000)
# x        = (c(0, 0.33, 0.66) - mean(data$X[,1])) / sd(data$X[,1])
x        = c(0, 0.5, 0.9)
f_yx     = matrix(0, nrow = length(y_grid), ncol = length(x)) # (1000 x 3)

## compute conditional density for each value of x
source("simulation.R")
# store the conditional density of (y | x) for each covariate
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



# initialize parameters needed to pass into initPriors() function
# TODO: figure out the appropriate value to initialize pi_d with
# C&S implementation seems to normalize the probabilities over all features (?)

alpha_0 = rep(1 / K, K);                         # pi_k params
m_0 = 0; xi_0 = 0.3;                             # beta params             
pi_d = 0.3;                                      # prior prob of inclusion
a_0 = 1; b_0 = 1;                                # tau params
VERBOSE = TRUE; tol = 1e-3;  max_iter = 1e4;     # algo/converge params

# initialize objects needed for cavi
# test initPriors() function
source("initPr.R")
source("initVP.R")
prior = initPriors(y, X, K, 
                   alpha_0,                    # pi   : concentration param
                   m_0, xi_0,                  # beta : mean, scale for slab
                   pi_d,                       # beta : prior pip
                   a_0, b_0)                   # tau  : shape, rate

theta = initVarParams_indep(y, X, N, D, K)

# perform each of the updates:
source("updateFunctions.R")
theta_p = theta


## rnkUpdate() : r_nk values changed -------------------------------------------
theta_p = rnkUpdate(prior, theta_p) # pass in updated theta: theta_p

head(theta_p$r_nk)

## wtUpdate() : alpha_k values changed -----------------------------------------
theta_p = wtUpdate(prior, theta_p) # pass in updated theta: theta_p

theta_p$alpha_k

## ssUpdate() : m_d, Q_d values updated ----------------------------------------

# before update values
theta_p$lambda_d
theta_p$m_d
theta_p$Q_d_inv

theta_p = ssUpdate(prior, theta_p) # pass in updated theta: theta_p 

# after update values
theta_p$lambda_d
theta_p$m_d
theta_p$Q_d_inv

## betak_update() : m_k, Sigma_k values updated --------------------------------

# before update values
theta_p$m_k
theta_p$Sigma_k

theta_p = betak_update(prior, theta_p) # pass in updated theta: theta_p

# after update values
theta_p$m_k
theta_p$Sigma_k


## precUpdate() : a_k, b_k values updated --------------------------------------

# before update values
theta_p$a_k
theta_p$b_k

theta_p = precUpdate(prior, theta_p) # pass in updated theta: theta_p

theta_p$a_k 
theta_p$b_k # these values are massive (>10k)




for (t in 1:200) {
    theta_p = rnkUpdate(prior, theta_p)
    theta_p = wtUpdate(prior, theta_p) 
    theta_p = ssUpdate(prior, theta_p) 
    theta_p = betak_update(prior, theta_p)
    theta_p = precUpdate(prior, theta_p)
}
















