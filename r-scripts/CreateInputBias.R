#' Create a dataframe of inputs used to evaluate GP metamodel of Bias
#'
#' The input to GP metamodel of bias, includes also time points and 
#' axial locations (for TC and DP).
#' These two input variables are kept constant during calibration, while the 
#' other variables related to boundary conditions are sampled by MCMC in each
#' iteration
#' 
#' @param xx Numeric vector. BC-related model parameters values for GP 
#'  metamodel of Bias
#' @param str_names Character vector. Names of the input parameters
#' @param time_pts Numeric vector. The time points in the bias formulation
#' @param ax_locs Numeric vector. Axial locations in the bias formulation
#' @param idx_time_length Numeric vector. The number of time points per axial
#'  location 
#' @return Dataframe. Expanded grid of BC-related model parameters column 
#'  binded with time points and axial location (if specified) grid
CreateInputBias <- function(xx, str_names, time_pts, 
                            ax_locs = NULL, idx_time_length = NULL)
{
    if (!is.null(ax_locs))
    {
        n_col <- length(str_names)
        xx_temp <- matrix(0, nrow = sum(idx_time_length), ncol = n_col)
        k <- 1
        for (i in 1:length(ax_locs))
        {
            for (j in 1:idx_time_length[i])
            {
                xx_temp[k, ] <- c(xx, ax_locs[i], time_pts[j])
                k <- k + 1
            }
        }
    } else 
    {
        xx_temp <- matrix(xx, nrow = 1)
        xx_temp <- cbind(xx_temp[rep(1:nrow(xx_temp), length(time_pts)),], 
                         matrix(time_pts, nrow = length(time_pts)))
    }

    # return dataframe
    xx_df <- data.frame(xx_temp)
    names(xx_df) <- str_names
    
    return(xx_df)
}
