# reconstruct_mean.R
# Sandboxing: reconstruct mean of Gaussian Process
install.packages("mvtnorm")
library(mvtnorm)
library(DiceKriging)


xx_test <- data.frame(x1 = 0.5, x2 = 0.5, x3 = 0.5, x4 = 0.5,
                      x5 = 0.5, x6 = 0.5, x7 = 0.5, x8 = 0.5,
                      x9 = 0.5, x10 = 0.5, x11 = 0.5, x12 = 0.5)
num_pc <- 3

# Read Data ----
trc_data <- readRDS("../postpro-gp-training/result-compiled/febaTrans216Ext-febaVars12Influential-lhs_120_12_1.Rds")
# GP Models
trc_gps_tc <- readRDS("../trace-gp/result-gp/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-tc-pca-gp-matern5_2.Rds")
# PCA Model
trc_pca_tc <- readRDS("../trace-gp/result-pc/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-tc-pca.Rds")
# Bias Model (discrepancy)
trc_dis_tc <- readRDS("../trace-bias/result-bias/tc/febaTrans216Ext-febaVars4BC-sobol_32_4-bias-tc.Rds")

time_steps <- unique(trc_dis_tc@X[,6])[-1]
ax_locs <- unique(trc_dis_tc@X[,5])
xx_test_bias <- data.frame(x1 = c(), x2 = c(), x3 = c(), x4 = c(), 
                           z = c(), t = c())

for (i in 1:length(ax_locs)) 
{
    n_ts <- length(idx_time_restricted[[i]][-1])
    for (j in 1:n_ts)
    {
        xx_test_bias <- rbind(xx_test_bias, data.frame(x1 = xx_test[1],
                                                       x2 = xx_test[2],
                                                       x3 = xx_test[3],
                                                       x4 = xx_test[4],
                                                        z = ax_locs[i],
                                                        t = time_steps[j]))    
    }
}

pred_mu <- c()
pred_sd <- c()
for (i in 1:num_pc)
{
    pred_mu <- c(pred_mu, 
                 predict(trc_gps_tc[[i]], newdata = xx_test, "SK")$mean)
    pred_sd <- c(pred_sd, 
                 predict(trc_gps_tc[[i]], newdata = xx_test, "SK")$sd)
}
pred_score <- matrix(pred_score, nrow = 3)
pred_var <- diag(pred_sd^2)

(trc_pca_tc$yy_train_loadings[,1:num_pc] * pred_score)[,3] == trc_pca_tc$yy_train_loadings[,num_pc] * pred_score[3]
mean_vector <- trc_pca_tc$yy_train_ave + trc_pca_tc$yy_train_loadings[,1:num_pc] %*% pred_score
var_matrix  <- trc_pca_tc$yy_train_loadings[,1:num_pc] %*% pred_var %*% t(trc_pca_tc$yy_train_loadings[,1:num_pc])

dim(trc_pca_tc$yy_train_loadings[,1:num_pc])
dim(pred_score)

bias_var_matrix <- covMatrix(trc_dis_tc@covariance, as.matrix(xx_test_bias))$C

length(time_steps)




# Create prototypical function ----
reconstruct_mean <- function(params, num_pc, trc_gps, trc_pca)
{
    
}

plot(trc_data$time[-1][idx_time_restricted_vector], matrix(mean_vector[idx_time_restricted_vector], ncol = 8)[,6])
lines(trc_data$time[-1], matrix(trc_pca_tc$yy_train_ave, ncol = 8)[,6])
lines(trc_data$time, trc_data$nominal[,6], col = "blue")


var_matrix <- trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,1:num_pc] %*% pred_var %*% t(trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,1:num_pc])
var_matrix_truncated <- trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,(num_pc+1):120] %*% diag(replicate(1,117)) %*% t(trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,(num_pc+1):120])

log(
    dmvnorm(tc_exp_vector, 
    mean = mean_vector[idx_time_restricted_vector]*1.1, 
    sigma = var_matrix + var_matrix_truncated + bias_var_matrix)
    )

# Create data vector
tc_exp_vector <- c()
for (j in 1:8)
{
    tc_exp_vector <- c(tc_exp_vector,
                       trc_data$exp_data[[1]][2:idx_exp_tquench[j],j+1])
}
tc_exp_vector <- tc_exp_vector + 273.15


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
plot(lala[,1], lala[,2])

