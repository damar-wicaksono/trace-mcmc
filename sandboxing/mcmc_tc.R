# mcmc_tc.R
# Sandboxing: Compute the likelihood of temperature prediction
install.packages("devtools")
library(devtools)
install_github("LaplacesDemonR/LaplacesDemon")
library(LaplacesDemon)

# 1. Read Data
# 1.a. Full TRACE Output of a given sample size, training design, replication
trc_data <- readRDS("../postpro-gp-training/result-compiled/febaTrans216Ext-febaVars12Influential-lhs_120_12_1.Rds")
# 1.b. GP Models of Standardized PC Scores
trc_gps_tc <- readRDS("../trace-gp/result-gp/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-tc-pca-gp-matern5_2.Rds")
# 1.c. PCA Result
trc_pca_tc <- readRDS("../trace-gp/result-pc/febaTrans216Ext-febaVars12Influential-lhs_120_12_1-tc-pca.Rds")
# 1.d. GP Model of Bias
trc_dis_tc <- readRDS("../trace-bias/result-bias/tc/febaTrans216Ext-febaVars4BC-sobol_32_4-bias-tc.Rds")

num_pc_tc <- 5

# 2. Get Restricted Time
# 2.a. TC Output (Up to Quenching)
# 2.b. CO Output (Up to Saturation)
# 3. Construct experimental data based on the restricted time
# 4. Construct data frame for sampled inputs, pc scores
# 5. Construct data frame for sampled inputs, bias
# 6. Compute the mean vector based on the restricted time
ave_vec_lmc <- function(xx_test, num_pc, idx_time, trc_gps, trc_pca)
{
    pred_mu <- c()
    for (i in 1:num_pc)
    {
        pred_mu <- c(pred_mu, 
                     predict(trc_gps[[i]], newdata = xx_test, "SK")$mean)
    }
    pred_mu <- matrix(pred_mu, nrow = num_pc)

    ave_vec <- trc_pca$yy_train_ave + trc_pca$yy_train_loadings[,1:num_pc] %*% pred_mu
    ave_vec <- ave_vec[idx_time]
    
    return(ave_vec)
}
ave_vec_tc <- ave_vec_lmc(xx_test, num_pc_tc, idx_time_restricted_vector, trc_gps_tc, trc_pca_tc)

# 7.a. Kriging Variance
# 7. Compute the variance matrix based on the restricted time
var_mat_lmc <- function(xx_test, num_pc, idx_time, trc_gps, trc_pca)
{
    n_pcs <- ncol(trc_pca$yy_train_loadings)
    # Kriging variance
    pred_sd <- c()
    for (i in 1:num_pc)
    {
        pred_sd <- c(pred_sd, 
                     predict(trc_gps[[i]], newdata = xx_test, "SK")$sd)
    }
    pred_var <- diag(pred_sd^2)
    
    var_mat_krig <- trc_pca$yy_train_loadings[idx_time, 1:num_pc] %*% 
        pred_var %*% t(trc_pca$yy_train_loadings[idx_time, 1:num_pc_tc])
    
    # PC Truncation Variance
    var_mat_trun <- trc_pca$yy_train_loadings[idx_time, (num_pc+1):n_pcs] %*% 
        diag(replicate(1, n_pcs-num_pc)) %*% 
        t(trc_pca$yy_train_loadings[idx_time, (num_pc+1):n_pcs])
    
    return(var_mat_krig + var_mat_trun)
    
}
diag(var_mat_lmc(xx_test, num_pc_tc, idx_time_restricted_vector, trc_gps_tc, trc_pca_tc))
diag(var_mat_krig_tc + var_mat_trun_tc)

var_mat_dis <- function(xx_test_bias, var_param, trc_dis_gp)
{
    # Create a new covariance 
    covtype <- trc_dis_gp@covariance@name
    d <- trc_dis_gp@covariance@d
    var_names <- trc_dis_gp@covariance@var.names
    if (covtype == "powexp")
    {
        coef_cov <- c(trc_dis_gp@covariance@range.val, 
                      trc_dis_gp@covariance@shape.val)    
    } else {
        coef_cov <- trc_dis_gp@covariance@range.val
    }
    if (is.na(var_param))
    {
        var_param <- trc_dis_gp@covariance@sd2
    } 
    
    cov_bias <- covStruct.create(covtype = covtype,
                                 d = d,
                                 known.covparam = "All",
                                 var.names = var_names,
                                 coef.cov = coef_cov,
                                 coef.var = var_param)
    # Discrepancy Variance
    var_mat_dis <- covMatrix(cov_bias, 
                             as.matrix(xx_test_bias))$C
    
    return(var_mat_dis)
}
diag(var_mat_dis(xx_test_bias_tc, NA, trc_dis_tc))^0.5
diag(covMatrix(trc_dis_tc@covariance, as.matrix(xx_test_bias_tc))$C)^0.5

# 7.b. PC Truncation Variance
# 7.c. Bias Variance
# 8. Given exp data vector, mean vector, and variance matrix compute likelihood
normal_ll <- function(exp_vec, ave_vec, var_mat)
{
    # Inverse by SVD
    svd_var_mat <- svd(var_mat)
    inv_var_mat <- svd_var_mat$v %*% diag(1/svd_var_mat$d) %*% t(svd_var_mat$u)
    
    k <- length(exp_vec)
    b1 <- matrix(exp_vec - ave_vec, ncol = k) %*% inv_var_mat %*% 
        matrix(exp_vec - ave_vec, nrow = k)
    b2 <- determinant(var_mat, logarithm = T)$modulus[1]
    b3 <- k * log(2*pi)
    
    return(-0.5*(b1 + b2 + b3))
}

var_mat <- var_mat_lmc(xx_test, num_pc_tc, idx_time_restricted_vector, trc_gps_tc, trc_pca_tc) + var_mat_dis(xx_test_bias_tc, NA, trc_dis_tc)
normal_ll(tc_exp_vector, ave_vec_tc, var_mat)

# Implementing Laplace Demon Model

# Set up model parameters
N <- nrow(tc_exp_vector)    # Number of data points
y <- tc_exp_vector          # Observed data
J <- ncol(xx_test)
param_names <- as.parm.names(list(x = rep(0,12), sigma = 0))
pos_x <- grep("x", param_names)
pos_sigma <- grep("sigma", param_names)
# Set up parameter generating function
# Creating list of MyData
MyFEBAData <- list(
    J = J,
    mon.names = "LP",
    parm.names = param_names,
    pos.x = pos_x,
    pos.sigma = pos_sigma,
    y = y
)

# Setup the Model, Multivariate Normal
Model <- function(params, Data)
{
    # Set up inputs
    x <- interval(params[Data$pos.x], 1e-100, 1)
    #x <- params[Data$pos_x]
    x_bias <- create_xx_bias(x[1:4], trc_dis_tc)
    params[Data$pos.x] <- x
    x <- create_xx_test(x)
    sigma <- interval(params[Data$pos.sigma], 1e-100, Inf)
    params[Data$pos.sigma] <- sigma
    
    # Log-Prior
    x_prior <- dunif(as.numeric(x), 0, 1, log = T)
    sigma_prior <- dhalfcauchy(sigma, 1000, log = T)

    # Log-Likelihood
    ave_vec <- ave_vec_lmc(x, num_pc_tc, idx_time_restricted_vector, 
                           trc_gps_tc, trc_pca_tc)
    var_mat <- var_mat_lmc(x, num_pc_tc, 
                           idx_time_restricted_vector, trc_gps_tc, trc_pca_tc) + 
        var_mat_dis(x_bias, sigma, trc_dis_tc)
    LL <- normal_ll(Data$y, ave_vec, var_mat)
    
    # Log-Posterior
    LP <- LL #+ sum(x_prior) + sigma_prior
    
    # Construct list
    Modelout <- list(
        LP = LP,
        Dev = -2 * LL,
        Monitor = LP,
        yhat = rmvnorm(1, ave_vec, var_mat),
        parm = params
    )
    return(Modelout)
}
Model(params=c(rep(0.5,12),2478), MyFEBAData)

Fit <- LaplacesDemon(Model, Data=MyFEBAData, c(rep(0.5,12), 1000),
                     Covar=NULL, Iterations=5000, Status=100, Thinning=1,
                     Algorithm="AFSS", Specs=list(A=Inf, B=NULL, m=100, n=0, w=1))




dunif(as.numeric(x_test), 0, 1, log = T)
x_test <- create_xx_test(rep(0,12))
x_test
var_mat_lmc(x_test, num_pc_tc, idx_time_restricted_vector, trc_gps_tc, trc_pca_tc)

create_xx_test <- function(x)
{
    xx_test <- data.frame(matrix(0, ncol = length(x), nrow = 0))
    x_names <- c()
    for (i in 1:length(x)) x_names <- c(x_names, paste("x", i, sep = ""))
    colnames(xx_test) <- x_names
    xx_test[1,] <- x
    return(xx_test)
}


create_xx_bias <- function(xx_bias, trc_dis_gp)
{
    time_steps <- unique(trc_dis_gp@X[,6])[-1]
    ax_locs <- unique(trc_dis_gp@X[,5])
    xx_bias_complete <- data.frame(x1 = c(), x2 = c(), x3 = c(), x4 = c(), 
                                   z = c(), t = c())
    for (i in 1:length(ax_locs)) 
    {
        n_ts <- length(idx_time_restricted[[i]][-1])
        for (j in 1:n_ts)
        {
            xx_bias_complete <- rbind(xx_bias_complete, 
                                      data.frame(x1 = xx_bias[1],
                                                 x2 = xx_bias[2],
                                                 x3 = xx_bias[3],
                                                 x4 = xx_bias[4],
                                                  z = ax_locs[i],
                                                  t = time_steps[j]))    
        }
    }
    return(xx_bias_complete)
}
create_xx_bias(xx_test[1:4], trc_dis_tc)
