

setwd(getwd())

source("globals.R")

setwd(HOME_DIR)
source(DP_BDR) 
source(DENSITY)


N = 500
K = 2
synth_data_1d = r_dpmix2(N)

y      = synth_data_1d$y
X      = as.matrix((synth_data_1d$X)[,2])
# X      = synth_data_1d$X
params = synth_data_1d$param_mat
y_grid = seq(-3, 1.5, length.out = 1000)

# 6th sublot in Figure 3
df_yx = data.frame(x = X[,1], y = y)
# df_yx = data.frame(x = X[,2], y = y)
ggplot(df_yx, aes(x = x, y = y)) + geom_point() + stat_smooth() + theme_bw()


# 1st subplot

x         = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
f_yx      = matrix(0, nrow = length(y_grid), ncol = length(x))

for(i in 1:length(x)) {
    # f_yx[,i] = d_dpmix2(y_grid, x[i], -1 + 2 * x[i])
    f_yx[,i] = d_dpmix2(y_grid = y_grid, x = c(x[i]), 
                        params = c(x[i], -2 * x[i]))
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
source(paste(COV_DEP, VARBDR, sep = '/'))    # load cov-dep varbdr.R
theta1_2 = varbdr(y, X, K)                   # store CAVI results


# save variational parameters
# saveRDS(theta1_1, file = "dp_results/theta_indep_dp2_K2_N1000.RDS")
saveRDS(theta1_2, file = "dp_results/theta_2.27_N_1e4_int.RDS")


# old1_1 = readRDS(file = "dp_results/theta_indep_dp2_K2_N1000.RDS")
# old1_2 = readRDS(file = "dp_results/theta_dep_dp2_K2_N1000.RDS")

plotELBO(theta1_2)





# generate plots ---------------------------------------------------------------

# theta1      = list(old1_2)      # list of var. params for each alg
theta1      = list(theta1_2)

# generate plot overlays of the percentiles
source(DENSITY)
x = c(0.15, 0.25, 0.49, 0.75, 0.88, 0.95)
approx_d  = list(py_bouch)            # list of approx density functions
den_label = c("cov-dep (b)")          # labels for each approx density
p_list2 = xQuantileDensity(x, params, d_dpmix2,
                           approx_d, den_label, theta1, K,
                           y_grid)

multiplot(plotlist = as.list(sapply(p_list2, function(x) x[1])), cols = 3)





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
