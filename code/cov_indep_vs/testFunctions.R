

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
X = scale(X, scale = TRUE)

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





