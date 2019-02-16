
setwd(getwd())

source("globals.R")

setwd(HOME_DIR)
source(UNI_CDE) 
source(DENSITY)


N = 500
K = 3

synth_data_bimodal = r_binorm(N)
y      = synth_data_bimodal$y
X      = synth_data_bimodal$X
params = synth_data_bimodal$param_mat
y_grid = seq(-5, 5, length.out = 1000)

x         = c(-0.6745, 0, 0.6745)
f_yx      = matrix(0, nrow = length(y_grid), ncol = length(x))

for(i in 1:length(x)) {
    f_yx[,i] = d_binorm(y_grid = y_grid, x = c(1, x[i]), 
                        params = c(x[i] - 1.5, x[i] + 1.5))
}

df_y      = data.frame(x = y_grid, f_yx)
df_y_long = melt(df_y, id.vars = "x")
levels(df_y_long$variable) = c("x = -0.6745", "x = 0", "x = 0.6745")
ggplot(df_y_long, aes(x, y = value)) + geom_point(size = 0.9) +
    scale_x_continuous(breaks = seq(min(y_grid), max(y_grid), by = 1)) + 
    facet_wrap(~variable)


# run cavi ---------------------------------------------------------------------

# (1.1) covariate-INDEPENENT vb
source(paste(COV_INDEP, VARBDR, sep = '/'))      # load cov-indep varbdr.R
theta1_1 = varbdr(y = y, X = X, K = 3)               # store CAVI results

# (1.2) covariate-DEPENDENT vb
source(paste(COV_DEP, VARBDR, sep = '/'))        # load cov-dep varbdr.R
theta1_2 = varbdr(y = y, X = X, K = 3)               # store CAVI results


theta2      = list(theta1_2)      # list of var. params for each alg
DATA_GEN_ID = BIMODAL1


# generate plots ---------------------------------------------------------------


# generate plot overlays of the percentiles1
source(DENSITY)
approx_d  = list(py_bouch)            # list of approx density functions
den_label = c("cov-dep (b)")   # labels for each approx density
# params argument isn't used, function will redefine what it needs
# params as is corresponds to the parameters for the randomly generated x
# here, we care about the parameters associated with the quantiles
x         = c(-0.6745, 0, 0.6745)
y_grid = seq(-5, 5, length.out = 1000)
p_list3 = xQuantileDensity(x, params, d_binorm,
                           approx_d, den_label, theta2, K,
                           DATA_GEN_ID, y_grid)

multiplot(plotlist = as.list(sapply(p_list3, function(x) x[1])), cols = 3)

# compare with np results ------------------------------------------------------
library("np")
X = data.frame(X[,2])
y = data.frame(y)
xy_df = data.frame(y = y, x = X)
names(xy_df) = c("x", "y")

fy.x = npcdens(y ~ x, xy_df)

overlays = vector("list", length(x))
y_grid = seq(-5, 5, length.out = 1000)

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



