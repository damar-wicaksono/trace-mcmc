# reconstruct_mean.R
# Sandboxing: reconstruct mean of Gaussian Process
install.packages("mvtnorm")
library(mvtnorm)
library(DiceKriging)

d_ts <- 0.1
xx_test <- data.frame(x1 = 0.5, x2 = 0.5, x3 = 0.5, x4 = 0.5,
                      x5 = 0.5, x6 = 0.5, x7 = 0.5, x8 = 0.5,
                      x9 = 0.5, x10 = 0.5, x11 = 0.5, x12 = 0.5)
num_pc_dp <- 4

# Read Data ----
trc_data <- readRDS("../postpro-gp-training/result-compiled/febaTrans216Ext-febaVars12Influential-lhs_120_12_1.Rds")
trc_data_dp <- readRDS("../postpro-gp-training/result-compiled/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-dp_smoothed.Rds")
# GP Models
trc_gps_dp <- readRDS("../trace-gp/result-gp/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-dp_smoothed-pca-gp-matern5_2.Rds")
# PCA Model
trc_pca_dp <- readRDS("../trace-gp/result-pc/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-dp_smoothed-pca.Rds")
# Bias Model (discrepancy)
trc_dis_dp <- readRDS("../trace-bias/result-bias/dp/febaTrans216Ext-febaVars4BC-sobol_32_4-bias-dp.Rds")

# Get restricted time
idx_time_exp_dp <- trc_data$exp_data[[2]][,1] / d_ts + 1
idx_time_trc_dp <- c()
for (j in 1:4)
{
    idx_time_trc_dp <- c(idx_time_trc_dp, idx_time_exp_dp + (j - 1) * 10000)
}

time_steps <- unique(trc_dis_dp@X[,6])[-1]
ax_locs <- unique(trc_dis_tc@X[,5])
xx_test_bias_dp <- data.frame(x1 = c(), x2 = c(), x3 = c(), x4 = c(), 
                           z = c(), t = c())

for (i in unique(trc_dis_dp@X[,5])) 
{
    for (j in unique(trc_dis_dp@X[,6]))
    {
        xx_test_bias_dp <- rbind(xx_test_bias_dp, data.frame(x1 = xx_test[1],
                                                             x2 = xx_test[2],
                                                             x3 = xx_test[3],
                                                             x4 = xx_test[4],
                                                              z = i,
                                                              t = j))    
    }
}
dim(xx_test_bias_dp)

pred_mu_dp <- c()
pred_sd_dp <- c()
for (i in 1:num_pc_dp)
{
    pred_mu_dp <- c(pred_mu_dp, 
                    predict(trc_gps_dp[[i]], newdata = xx_test, "SK")$mean)
    pred_sd_dp <- c(pred_sd_dp, 
                    predict(trc_gps_dp[[i]], newdata = xx_test, "SK")$sd)
}
pred_mu_dp  <- matrix(pred_mu_dp, nrow = num_pc_dp)
pred_var_dp <- diag(pred_sd_dp^2)

ave_vec_dp <- trc_pca_dp$yy_train_ave  + trc_pca_dp$yy_train_loadings[,1:num_pc_dp] %*% pred_mu_dp
ave_vec_dp <- ave_vec_dp[idx_time_trc_dp]
var_mat_kriging_dp <- trc_pca_dp$yy_train_loadings[idx_time_trc_dp,1:num_pc_dp] %*% pred_var_dp %*% t(trc_pca_dp$yy_train_loadings[idx_time_trc_dp,1:num_pc_dp])
var_mat_truncat_dp <- trc_pca_dp$yy_train_loadings[idx_time_trc_dp,(num_pc_dp+1):120] %*% diag(replicate(1,116)) %*% t(trc_pca_dp$yy_train_loadings[idx_time_trc_dp,(num_pc_dp+1):120])
var_mat_discrep_dp <- covMatrix(trc_dis_dp@covariance, 
                                as.matrix(xx_test_bias_dp))$C

diag(var_mat_kriging_dp)^0.5
diag(var_mat_truncat_dp)^0.5
diag(var_mat_discrep_dp)^0.5

# Create data vector
dp_exp_vector <- c()
for (j in 1:4)
{
    dp_exp_vector <- c(dp_exp_vector,
                       trc_data$exp_data[[2]][,j+1])
}
dp_exp_vector <- dp_exp_vector * 1e5

log(
    dmvnorm(dp_exp_vector, 
            mean = ave_vec_dp, 
            sigma = var_mat_kriging_dp + var_mat_truncat_dp + var_mat_discrep_dp)
)

# Create prototypical function ----
reconstruct_mean <- function(params, num_pc, trc_gps, trc_pca)
{
    
}

plot(trc_data$time, trc_data_dp$nominal[,2], col = "blue", type = "l")
points(trc_data$exp_data[[2]][,1], 
       trc_data$exp_data[[2]][,3]*1e5, col = "green", pch = 2, lwd = 2)
points(trc_data$time[idx_time_exp_dp], 
       matrix(ave_vec_dp, ncol = 4)[,2], lwd = 2, pch = 4, col = "black")

lines(trc_data$time, matrix(trc_pca_dp$yy_train_ave, ncol = 4)[,4])


plot(trc_data$exp_data[[2]][,1], trc_data$exp_data[[2]][,3]*1e5, col = "blue")
for (i in 1:120) points(trc_data$time[idx_time_exp_dp], 
                        trc_data_dp$replicates[idx_time_exp_dp,i,4], col = "blue")


# Things to do:
# 1. Read Data
# 1.a. Full TRACE Output of a given sample size, training design, replication
# 1.b. GP Models of Standardized PC Scores
# 1.c. PCA Result
# 1.d. GP Model of Bias
# 2. Get Restricted Time
# 2.a. TC Output (Up to Quenching)
# 2.b. CO Output (Up to Saturation)
# 3. Construct experimental data based on the restricted time
# 4. Construct data frame for sampled inputs, pc scores
# 5. Construct data frame for sampled inputs, bias
# 6. Compute the mean vector based on the restricted time
# 7. Compute the variance matrix based on the restricted time
# 7.a. Kriging Variance
# 7.b. PC Truncation Variance
# 7.c. Bias Variance
# 8. Given exp data vector, mean vector, and variance matrix compute likelihood

plot(seq(-3,3, length.out = 100),
     dnorm(seq(-3,3, length.out = 100), mean = 0, sd = 1))
dv
lala <- rmvnorm(n=500, mean=c(1,2), sigma=matrix(c(4,2,2,3), ncol=2), method="chol")
plot(trc_gps_dp[[5]])

