# ll_tc.R
# Sandboxing: Compute the likelihood of temperature prediction
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

# Compute Restricted Time


# Set up input
d_ts <- 0.1
xx_test <- data.frame(x1 = 0.5, x2 = 0.5, x3 = 0.5, x4 = 0.5,
                      x5 = 0.5, x6 = 0.5, x7 = 0.5, x8 = 0.5,
                      x9 = 0.5, x10 = 0.5, x11 = 0.5, x12 = 0.5)
rep <- 101
xx_test <- xx_train[rep,]
xx_name <- c()
for (i in 1:12) xx_name <- c(xx_name, paste("x", i, sep = ""))
names(xx_test) <- xx_name

num_pc_tc <- 5

time_steps <- unique(trc_dis_tc@X[,6])[-1]
ax_locs <- unique(trc_dis_tc@X[,5])
xx_test_bias_tc <- data.frame(x1 = c(), x2 = c(), x3 = c(), x4 = c(), 
                              z = c(), t = c())

for (i in 1:length(ax_locs)) 
{
    n_ts <- length(idx_time_restricted[[i]][-1])
    for (j in 1:n_ts)
    {
        xx_test_bias_tc <- rbind(xx_test_bias_tc, data.frame(x1 = xx_test[1],
                                                             x2 = xx_test[2],
                                                             x3 = xx_test[3],
                                                             x4 = xx_test[4],
                                                              z = ax_locs[i],
                                                              t = time_steps[j]))    
    }
}

pred_mu_tc <- c()
pred_sd_tc <- c()
for (i in 1:num_pc_tc)
{
    pred_mu_tc <- c(pred_mu_tc, 
                    predict(trc_gps_tc[[i]], newdata = xx_test, "SK")$mean)
    pred_sd_tc <- c(pred_sd_tc, 
                    predict(trc_gps_tc[[i]], newdata = xx_test, "SK")$sd)
}
pred_mu_tc <- matrix(pred_mu_tc, nrow = num_pc_tc)
pred_var_tc <- diag(pred_sd_tc^2)

(trc_pca_tc$yy_train_loadings[,1:num_pc] * pred_score)[,3] == trc_pca_tc$yy_train_loadings[,num_pc] * pred_score[3]

ave_vec_tc <- trc_pca_tc$yy_train_ave + trc_pca_tc$yy_train_loadings[,1:num_pc_tc] %*% pred_mu_tc
ave_vec_tc <- ave_vec_tc[idx_time_restricted_vector]
var_mat_krig_tc <- trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,1:num_pc_tc] %*% 
    pred_var_tc %*% t(trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,1:num_pc_tc])
var_mat_trun_tc <- trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,(num_pc_tc+1):120] %*% 
    diag(replicate(1,120-num_pc_tc)) %*% t(trc_pca_tc$yy_train_loadings[idx_time_restricted_vector,(num_pc_tc+1):120])
var_mat_disc_tc <- covMatrix(trc_dis_tc@covariance, 
                             as.matrix(xx_test_bias_tc))$C

diag(var_mat_krig_tc)^0.5
diag(var_mat_trun_tc)^0.5
diag(var_mat_disc_tc)^0.5

# Create data vector
tc_exp_vector <- c()
for (j in 1:8)
{
    tc_exp_vector <- c(tc_exp_vector,
                       trc_data$exp_data[[1]][2:idx_exp_tquench[j],j+1])
}
tc_exp_vector <- tc_exp_vector + 273.15

# Log likelihood
log(
    dmvnorm(tc_exp_vector, 
            mean = ave_vec_tc, 
            sigma = var_mat_krig_tc + var_mat_trun_tc + var_mat_disc_tc)
)

# Slice output for a given axial location
j <- 5
slice <- 0
for (i in 1:(j-1))
{
    slice <- slice + length(idx_time_restricted[[i]][-1])
}
slice <- c((slice + 1):(slice + length(idx_time_restricted[[j]][-1])))

plot(trc_data$exp_data[[1]][2:idx_exp_tquench[j],1], 
     trc_data$exp_data[[1]][2:idx_exp_tquench[j],j+1]+273.15, 
     col = "green", pch = 2, lwd = 2, ylim = c(300, 1400))
points(trc_data$time[idx_time_restricted[[j]][-1]], 
       trc_data$replicates[idx_time_restricted[[j]][-1],rep,j], 
       col = "blue")
points(trc_data$time[idx_time_restricted[[j]][-1]], 
       ave_vec_tc[slice], 
       lwd = 2, pch = 4, col = "black")
points(trc_data$time[idx_time_restricted[[j]][-1]], 
       ave_vec_tc[slice] + 2 * ((diag(var_mat_krig_tc)^0.5)[slice] + (diag(var_mat_trun_tc)^0.5)[slice] + (diag(var_mat_disc_tc)^0.5)[slice]), 
       lwd = 2, pch = 6, col = "black")
points(trc_data$time[idx_time_restricted[[j]][-1]], 
       ave_vec_tc[slice] - 2 * ((diag(var_mat_krig_tc)^0.5)[slice] + (diag(var_mat_trun_tc)^0.5)[slice] + (diag(var_mat_disc_tc)^0.5)[slice]), 
       lwd = 2, pch = 6, col = "black")


lines(trc_data$time[idx_time_restricted[[8]][-1]], 
      matrix(mean_vector[idx_time_restricted_vector], ncol = 8)[,6])
lines(trc_data$time[-1], matrix(trc_pca_tc$yy_train_ave, ncol = 8)[,6])
lines(trc_data$time, trc_data$nominal[,6], col = "blue")


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

