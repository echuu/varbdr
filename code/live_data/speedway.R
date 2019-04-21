


setwd("C:/Users/chuu/varbdr/code")

source("globals.R")

setwd(HOME_DIR)
source(UNI_CDE) 
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation


setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")                              # load debug functions
sourceCpp("getVarParams.cpp")                        # load C++ implementation
source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation



library(hdrcde)
library(np)

# 
# lane2 : (1318 x 2) lane 2 of the 4-lane California freeway I-880
# flow  : traffic flow in vehicles per lane per hour
# speed : speed in miles per hour

plot(lane2) # plot of traffic speed vs. traffic flow

lane2_sub = lane2 %>% filter(flow >= 1000, flow <= 1620)

X = cbind(as.matrix(lane2$flow))
X = X / 100
X = cbind(1, X)
y = lane2$speed
N = nrow(X)
D = ncol(X)
K = 2
intercept = TRUE
max_iter  = 5e4

y_grid = seq(min(y), max(y), length.out = 2000)
x_in = c(1400)

meanParams = generateParams(y, X, N, D, K, intercept = FALSE, max_iter)
(m_k  = meanParams$m_k)
(mu_k = meanParams$mu_k)
theta_R  = varbdr(y = y, X = X, K, intercept = FALSE, m_k = m_k, mu_k = mu_k)    

theta_sf = varbdr_cpp(y, X, N, D, K, intercept, max_iter)

N_evals    = length(y_grid)                          # number of evaluations
data_ygrid = list(y = y_grid, x = c(1, x_in / 100))  # fixed for all evaluations

approx_density = py_bouch(theta_sf, K, data_ygrid) 

approx_df = data.frame(y = y_grid, p_y = approx_density)

ggplot(approx_df, aes(x = y, y = p_y)) + 
    geom_line(size = 0.8) + labs(x = "y", y = "p(y)") + theme_bw() 















###### debugging

m_0 = c(colMeans(X))                         # normal params              
Lambda_0 = diag(rep(1, ncol(X))) 
a_0 = 1
b_0 = 1                             # gamma params
g_0 = 0 
Sigma_0 = diag(rep(1, ncol(X)))


###### debugging

source(paste(COV_DEP, VARBDR, sep = '/'))            # load R implementation

prior = initPriors(y, X, K, m_0, Lambda_0, a_0, b_0, g_0, Sigma_0,
                   max_iter, tol, VERBOSE)

theta = initVarParams(y, X, N, D, K, intercept, max_iter, m_k, mu_k)

theta = eStep(theta, prior)
theta = mStep(theta, prior)
theta = elbo(theta, prior)

# check xi (previously had too large values), phi (overflow in exponential term)
head(theta$xi, 50)
head(theta$phi, 50)

# check the expectations
theta$e_ln_p_y
theta$e_ln_p_z
theta$e_ln_p_gamma
theta$e_ln_p_beta_tau
theta$e_ln_q_z
theta$e_ln_q_beta_tau
theta$e_ln_q_gamma



