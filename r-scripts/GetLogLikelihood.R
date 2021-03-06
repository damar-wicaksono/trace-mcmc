#' Wrapper function to compute the log-likelihood of a given parameters values
#'
#' @param x Numeric vector. Sampled model parameter, 13 of them. The last one
#'  is always the scale parameter of the bias model.
#' @param sampled_idx Numeric vector. Model parameters to be sampled by MCMC.
#' @param fix_bc Logical. Flag whether to fix the BC-related parameters in the
#'  GP metamodel
#' @param correct_bias Logical. Flag whether to correct for the bias
#' @param fix_scale Logical. Flag whether to fix the scale parameter of the 
#'  bias GP model
#' @param exp_vec Numeric vector. The vector of experimental data.
#' @param num_pc Numeric. The number of principal components for reconstructing
#'  full model output.
#' @param time_idx Numeric vector. The subset of the full model output to be
#'  compared with experimental data.
#' @param trc_gps_pcs List. List of km objects, the metamodels of PC scores.
#' @param trc_pca_ave Numeric vector. The average of full model output.
#' @param trc_pca_lds Matrix. The principal component loadings.
#' @param trc_gp_bias \code{km} objects. The GP model of bias.
#' @param ax_locs Numeric vector. Unique number representing axial locations.
#' @param time_pts Numeric vector. Unique number of time_pts.
#' @param exp_idx_max Numeric vector. The maximum number of experimental 
#'  data points per axial location.
GetLogLikelihood <- function(x, sampled_idx, fix_bc, correct_bias, fix_scale,
                             exp_vec, num_pc, time_idx,
                             trc_gps_pcs, trc_pca_ave, trc_pca_lds, 
                             trc_gp_bias, ax_locs, time_pts, exp_idx_max)
{
    xx <- rep(0.5, 12)
    if (fix_scale)
    {
        xx[sampled_idx] <- x
    } else
    {
        xx[sampled_idx] <- x[1:(length(x)-1)]
    }
    
    # Construct data frame for sampled inputs, pc scores
    str_names <- trc_gps_pcs[[1]]@covariance@var.names
    xx_pcs <- CreateInputPCScores(xx[1:12], str_names)
    # Fix BC parameter in the GP metamodel
    if (fix_bc)
    {
        xx_pcs[1:4] <- rep(0.5, 4)
    }

    # Construct data frame for sampled inputs, bias
    str_names <- trc_gp_bias@covariance@var.names
    xx_bias <- CreateInputBias(xx[1:4], 
                              str_names, time_pts, 
                              ax_locs, exp_idx_max)

    # Compute PC Scores
    pc_scores <- CalcPCScores(xx_pcs, trc_gps_pcs, num_pc)

    # Compute the mean vector based on the restricted time
    ave_vec <- CalcAveVec(pc_scores$mu, 
                          trc_pca_ave, trc_pca_lds, time_idx)
    
    # Compute the mean of bias model
    if (correct_bias)
    {
        bias_vec <- predict(trc_gp_bias, newdata = xx_bias, "SK")$mean
        ave_vec <- ave_vec + bias_vec
    }
    
    # Compute the variance matrix due to PC scores prediction error
    kriging_var_mat <- CalcKrigingVarMat(pc_scores$sd, 
                                         trc_pca_lds, time_idx)

    # Compute the variance matrix due to PC truncation
    truncation_var_mat <- CalcTruncationVarMat(trc_pca_lds, 
                                               length(pc_scores$sd), 
                                               time_idx)
    
    # Compute the variance matrix due to model bias
    if (!fix_scale)
    {
        bias_var_mat <- CalcBiasVarMat(xx_bias, 
            x[length(x)]^2, 
            trc_gp_bias@covariance)
    } else 
    {
        bias_var_mat <- CalcBiasVarMat(xx_bias, 
        trc_gp_bias@covariance@sd2, 
        trc_gp_bias@covariance)
    }
    
    # Compute the variance matrix based on the restricted time
    var_mat <- kriging_var_mat + truncation_var_mat + bias_var_mat

    # Given exp data vector, mean vector, and variance matrix comp. likelihood
    ll <- CalcLogLikelihood(exp_vec, ave_vec, var_mat)
    
    return(as.numeric(ll))
}
