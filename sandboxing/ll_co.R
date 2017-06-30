# ll_co.R
# Sandboxing: Compute the likelihood of liquid carryover
install.packages("mvtnorm")
library(mvtnorm)
library(DiceKriging)

d_ts <- 0.1
xx_test <- data.frame(x1 = 0.5, x2 = 0.5, x3 = 0.5, x4 = 0.5,
                      x5 = 0.5, x6 = 0.5, x7 = 0.5, x8 = 0.5,
                      x9 = 0.5, x10 = 0.5, x11 = 0.5, x12 = 0.5)
rep <- 120
xx_test <- xx_train[rep,]
xx_name <- c()
for (i in 1:12) xx_name <- c(xx_name, paste("x", i, sep = ""))
names(xx_test) <- xx_name

num_pc_co <- 3

# Read Data ----
trc_data <- readRDS("../postpro-gp-training/result-compiled/febaTrans216Ext-febaVars12Influential-lhs_120_12_1.Rds")
# GP Models
trc_gps_co <- readRDS("../trace-gp/result-gp/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-co-pca-gp-matern5_2.Rds")
# PCA Model
trc_pca_co <- readRDS("../trace-gp/result-pc/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-co-pca.Rds")
# Bias Model (discrepancy)
trc_dis_co <- readRDS("../trace-bias/result-bias/co/febaTrans216Ext-febaVars4BC-sobol_32_4-bias-co.Rds")
# Training inputs
xx_train <- read.csv("../../wd41-thesis/simulation/gp-training/dmfiles/lhs_120_12_1.csv", header = F)

# Get restricted time
# Make local variables
idx_time_exp_co <- 1:(sum(trc_data$exp_data[[3]][,2] < 10.0) + 1) # Last time point before saturation
idx_time_trc_co <- trc_data$exp_data[[3]][idx_time_exp_co,1] / d_ts + 1

time_steps <- unique(trc_dis_co@X[,5])

xx_test_bias_co <- data.frame(x1 = c(), x2 = c(), x3 = c(), x4 = c(), t = c())

for (i in unique(trc_dis_co@X[,5])) 
{
    xx_test_bias_co <- rbind(xx_test_bias_co, data.frame(x1 = xx_test[1],
                                                         x2 = xx_test[2],
                                                         x3 = xx_test[3],
                                                         x4 = xx_test[4],
                                                          t = i))    
}

dim(xx_test_bias_co)

pred_mu_co <- c()
pred_sd_co <- c()
for (i in 1:num_pc_co)
{
    pred_mu_co <- c(pred_mu_co, 
                    predict(trc_gps_co[[i]], newdata = xx_test, "SK")$mean)
    pred_sd_co <- c(pred_sd_co, 
                    predict(trc_gps_co[[i]], newdata = xx_test, "SK")$sd)
}
pred_mu_co  <- matrix(pred_mu_co, nrow = num_pc_co)
pred_var_co <- diag(pred_sd_co^2)
pred_dis_co <- predict(trc_dis_co, newdata = xx_test_bias_co, "SK")$mean

ave_vec_co <- trc_pca_co$yy_train_ave  + trc_pca_co$yy_train_loadings[,1:num_pc_co] %*% pred_mu_co
ave_vec_co <- ave_vec_co[idx_time_trc_co]
var_mat_kriging_co <- trc_pca_co$yy_train_loadings[idx_time_trc_co,1:num_pc_co] %*% pred_var_co %*% t(trc_pca_co$yy_train_loadings[idx_time_trc_co,1:num_pc_co])
var_mat_truncat_co <- trc_pca_co$yy_train_loadings[idx_time_trc_co,(num_pc_co+1):120] %*% diag(replicate(1,117)) %*% t(trc_pca_co$yy_train_loadings[idx_time_trc_co,(num_pc_co+1):120])
var_mat_discrep_co <- covMatrix(trc_dis_co@covariance, 
                                as.matrix(xx_test_bias_co))$C
diag(predict(trc_dis_co, newdata = xx_test_bias_co, "SK", cov=T)$cov)^0.5
diag(var_mat_kriging_co)^0.5
diag(var_mat_truncat_co)^0.5
diag(var_mat_discrep_co)^0.5

# Create data vector
co_exp_vector <- trc_data$exp_data[[3]][idx_time_exp_co,2]

# Log Likelihood
log(
    dmvnorm(co_exp_vector, 
            mean = ave_vec_co, 
            sigma = var_mat_kriging_co + var_mat_truncat_co + var_mat_discrep_co)
)

# Create prototypical function ----
reconstruct_mean <- function(params, num_pc, trc_gps, trc_pca)
{
    
}


plot(trc_data$exp_data[[3]][idx_time_exp_co,1], 
     trc_data$exp_data[[3]][idx_time_exp_co,2], 
     col = "green", pch = 2, lwd = 2)
lines(trc_data$time[idx_time_trc_co], trc_data$replicates[idx_time_trc_co,rep,13], 
      col = "blue", type = "l")
points(trc_data$time[idx_time_trc_co], ave_vec_co, 
       lwd = 2, pch = 4, col = "black")
points(trc_data$time[idx_time_trc_co], 
       ave_vec_co + 2 * (diag(var_mat_kriging_co)^0.5 + diag(var_mat_truncat_co)^0.5), 
       lwd = 2, pch = 6, col = "black")
points(trc_data$time[idx_time_trc_co], 
       ave_vec_co - 2 * (diag(var_mat_kriging_co)^0.5 + diag(var_mat_truncat_co)^0.5), 
       lwd = 2, pch = 6, col = "black")
points(trc_data$time[idx_time_trc_co], ave_vec_co + pred_dis_co, 
       lwd = 2, pch = 6, col = "blue")
lines(trc_data$time, trc_data$replicates[,rep,13], 
      col = "black", type = "l")
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
# 6.a. Compute the mean from Kriging prediction of PC scores
# 6.b. Compute the mean bias if requested
# 7. Compute the variance matrix based on the restricted time
# 7.a. Kriging Variance
# 7.b. PC Truncation Variance
# 7.c. Bias Variance
# 8. Given exp data vector, mean vector, and variance matrix compute likelihood

