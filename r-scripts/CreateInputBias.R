#' Create a dataframe of inputs used to evaluate GP metamodel of Bias
#'
#' The input to GP metamodel of bias, includes also time points and 
#' axial locations (for TC and DP).
#' These two input variables are kept constant during calibration, while the 
#' other variables related to boundary conditions are sampled by MCMC. 
#' 
#' @param xx Numeric vector. BC-related model parameters values for GP 
#'  metamodel of Bias
#' @param aux_variables Matrix. BC 
#' @param str_name Character vector. Names of the input parameters
#' @return Dataframe. A single element dataframe to be used in GP metamodel of
#'  PC scores
CreateInputBias <- function(xx, str_names, time_points, ax_locs)
{
    if (!is.na(ax_locs))
    {
        length_xx <- unique(ax_locs) * unique(time_points)
        xx_df <- data.frame(matrix(0, ncol = length(xx) + 2, nrow = length_xx))
    } else 
    {
        length_xx <- unique(time_points)
        xx_df <- data.frame(matrix(0, ncol = length(xx) + 1, nrow = length_xx))
    }
    
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