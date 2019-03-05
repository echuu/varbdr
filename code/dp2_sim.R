

setwd(getwd())

source("globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)

N = 5e4
K = 2
synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = synth_data_1d$X
params = synth_data_1d$param_mat
y_grid = seq(-3, 1.5, length.out = 1000)

# 6th sublot in Figure 3
df_yx = data.frame(x = X[,1], y = y)
ggplot(df_yx, aes(x = x, y = y)) + geom_point() + stat_smooth() + theme_bw()


# 1st subplot
x         = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
f_yx      = matrix(0, nrow = length(y_grid), ncol = length(x))

for(i in 1:length(x)) {
    # f_yx[,i] = d_dpmix2(y_grid, x[i], -1 + 2 * x[i])
    f_yx[,i] = d_dpmix2(y_grid = y_grid, x = c(x[i]))
}

df_y      = data.frame(x = y_grid, f_yx)
df_y_long = melt(df_y, id.vars = "x")
levels(df_y_long$variable) = c("x = 0.15", "x = 0.25", "x = 0.49",
                               "x = 0.75", "x = 0.88", "x = 0.95")
ggplot(df_y_long, aes(x, y = value)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(min(y_grid), max(y_grid), by = 1)) + 
    facet_wrap(~variable)

# run cavi ---------------------------------------------------------------------

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))         # load cov-dep varbdr.R
source(paste(COV_DEP, FAST_VARBDR, sep = '/'))    # load cov-dep fast varbdr.R

start_time <- Sys.time()
theta0 = varbdr(y, X, K)                   # store CAVI results
end_time <- Sys.time()
(diff_fix = end_time - start_time)


microbenchmark(varbdr(y, X, K),
               fast_varbdr(y, X, K))


theta1 = fast_varbdr(y, X, K)


# save variational parameters
saveRDS(theta0, file = "dp_results/theta_N_500_K_2.RDS")

theta_500_2 = readRDS(file = "dp_results/theta_N_500_K_2.RDS") # N = 1e3, K = 2
theta_1e3_2 = readRDS(file = "dp_results/theta_N_1e3_K_2.RDS") # N = 1e3, K = 2
theta_1e3_3 = readRDS(file = "dp_results/theta_N_1e3_K_3.RDS") # N = 1e3, K = 3
theta_1e3_4 = readRDS(file = "dp_results/theta_N_1e3_K_4.RDS") # N = 1e3, K = 4
theta_1e4_2 = readRDS(file = "dp_results/theta_N_1e4_K_2.RDS") # N = 1e4, K = 2
theta_1e4_3 = readRDS(file = "dp_results/theta_N_1e4_K_3.RDS") # N = 1e4, K = 3
theta_2e4_2 = readRDS(file = "dp_results/theta_N_2e4_K_2.RDS") # N = 2e4, K = 2


# generate plots ---------------------------------------------------------------


overlayPlots(theta_1e3_4, 4)


overlayPlots = function(theta, K, den_label,
                        x = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)) {
    
    source(DENSITY)

    n = length(theta)
    approx_d  = rep(list(py_bouch), n)
    # approx_d = list(py_bouch, py_bouch, py_bouch, py_bouch)
    # den_label = c("cov-dep-1", "cov-dep-2", "cov-dep-3", "cov-dep-4")
    p_list2 = xQuantileDensity(x, params, d_dpmix2,
                               approx_d, den_label, theta, K,
                               y_grid)
    
    p = multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), 
                  cols = 3)
    return(p)
}


overlayPlots(theta = list(theta_500_2, theta_1e3_2, theta_1e4_2, theta_2e4_2), 
             K = 2, den_label = c("N=500", "N=1e3", 
                                  "N=1e4", "N=2e4"))

overlayPlots(list(theta_1e3_3, theta_1e4_3), K = 3,
             den_label = c("N=1e3", "N=1e4"))

overlayPlots(list(theta_1e3_4), K = 4, 
             den_label = c("N=1e4"))


# theta1      = list(old1_2)      # list of var. params for each alg
theta1_list    = list(theta0)

# generate plot overlays of the percentiles
source(DENSITY)
x = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
approx_d  = list(py_bouch, py_bouch)    # list of approx density functions
den_label = c("cov-dep")                # labels for each approx density
p_list2 = xQuantileDensity(x, params, d_dpmix2,
                           approx_d, den_label, theta1_list, K,
                           y_grid)

multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), cols = 3)

# -----------------------------------------------------------------------------










# compare with np results ------------------------------------------------------
library("np")
X = data.frame(X[,2])
y = data.frame(y)
xy_df = data.frame(y = y, x = X)
names(xy_df) = c("x", "y")

fy.x = npcdens(y ~ x, xy_df)

overlays = vector("list", length(x))
y_grid = seq(-0.5, 1.5, length.out = 1000)

for (i in 1:length(x)) {
    y_eval = data.frame(y = y_grid, x = x[i])
    fy_eval = predict(fy.x, newdata = y_eval)
    
    fy_df = data.frame(y = y_grid, np = fy_eval)
    
    # ggplot(fy_df, aes(x = y, y = np)) + geom_line()
    
    # overlay np estimate with the other approximations ------------------------
    
    # x = 0.15
    full_approx_df = merge(p_list2[[i]]$approx_df, fy_df, by = "y")
    approx_df_long = melt(full_approx_df, 
                          measure.vars = c(den_label, "true", "np"))
    
    overlays[[i]] = ggplot(approx_df_long, 
                           aes(x = y, y = value, colour = variable)) + 
        geom_line(size = 0.8) + labs(x = "y", y = "p(y)") + theme_bw() +
        theme(legend.position = "none")
}

multiplot(plotlist = overlays, cols = 3)
