
# make sure to use the correct constructor for initiVarParams()
source("C:/Users/chuu/varbdr/code/globals.R")

setwd(HOME_DIR)
source("plotCD.R")
source(DP_BDR) 
source(DENSITY)
source(paste(COV_DEP, VARBDR, sep = '/')) # load cov-dep varbdr.R

setwd("C:/Users/chuu/varbdr/code/cpp_code")
source("debug_funcs.R")
sourceCpp("getVarParams.cpp")

N = 1e3
D = 1
K = 2
synth_data_1d = r_dpmix2(N)
y = synth_data_1d$y
X = synth_data_1d$X
y_grid = seq(-3, 1.5, length.out = 1000)
intercept = FALSE
max_iter = 1e4

# plot y ~ X
df_yx = data.frame(x = X[,1], y = y)
yx_plot = ggplot(df_yx, aes(x = x, y = y)) + 
    geom_point(size = 1.1) + theme_bw()


# meanParams = generateParams(y, X, N, D, K, intercept, max_iter)
# m_k  = meanParams$m_k
# mu_k = meanParams$mu_k 


# run both cavi algorithms
# theta0_R = varbdr(y, X, K, m_k = m_k, mu_k = mu_k)   # store CAVI results

theta0_cpp = varbdr_cpp(y, X, N, D, K, intercept, max_iter)

# saveRDS(theta0_cpp, file = "dp_results/theta_N_1e5_K_2.RDS")
# check equality after convergence
# checkCorrectness(theta0_R, theta0_cpp)

# compare speed
# microbenchmark(varbdr(y, X, K, m_k = m_k, mu_k = mu_k),
#               varbdr_cpp(y, X, N, D, K, intercept, max_iter),
#               times = 25)


## check plots -- need to re-edit the functions to take C++ return
setwd(HOME_DIR)
theta = readRDS(file = "dp_results/theta_N_1e4_K_2.RDS")

theta = list(theta0_R)
overlayPlots = function(theta, K, den_label,
                        x = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)) {
    
    # source(DENSITY)
    
    n = length(theta)
    approx_d  = rep(list(py_bouch), n)
    p_list2 = xQuantileDensity(x, params, d_dpmix2,
                               approx_d, den_label, theta, K,
                               y_grid)
    
    p = multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), 
                  cols = 3)
    return(p)
}

overlayPlots(list(theta0_cpp), K = 2, den_label = c("N=500"))


## using updated plotting function
X = data.frame(X[,1])
y = data.frame(y)
xy_df = data.frame(y = y, x = X)
names(xy_df) = c("y", "x")

fy.x_dp = npcdens(y ~ x, xy_df)


x = c(0.03, 0.32, 0.49, 0.75, 0.97)
xlabs = paste("x = ", x, sep = "")
y_grid = seq(-3, 1.5, length.out = 1000)
dp2_cd = plotCD(theta, K, x, y_grid, xlabs, 
                true_den = d_dpmix2, k_den = fy.x_dp)
# multiplot(plotlist = dp2_cd$cd_plots, cols = 3)

## create visual for poster
multiplot(yx_plot, dp2_cd$cd_plots[[1]], dp2_cd$cd_plots[[5]],
          cols = 3)



