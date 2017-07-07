#' Create a dataframe of inputs used to evaluate GP metamodel of Bias
#'
#' The input to GP metamodel of bias, includes also time points and 
#' axial locations (for TC and DP).
#' These two input variables are kept constant during calibration, while the 
#' other variables related to boundary conditions are sampled by MCMC. 
#' 
#' @param xx Numeric vector. BC-related model parameters values for GP 
#'  metamodel of Bias
#' @param str_names Character vector. Names of the input parameters
#' @param time_pts Numeric vector.
#' @param ax_locs Numeric vector. BC 
#' @return Dataframe. Expanded grid of BC-related model parameters column 
#'  binded with time points and axial location (if specified) grid
CreateInputBias <- function(xx, str_names, time_pts, ax_locs = NULL)
{
    xx <- matrix(xx, nrow = 1)
    if (!is.null(ax_locs))
    {
        # Create Z-T grid but switch column because time is varying fastest
        zt_grid <- expand.grid(time_pts, ax_locs)[,c(2,1)] 
    } else 
    {
        zt_grid <- expand.grid(time_pts)
    }
    
    # Repeat xx as many as there is row in zt-grid before column binding
    xx_df <- cbind(xx[rep(1:nrow(xx), times = nrow(zt_grid)),], zt_grid)
    names(xx_df) <- str_names
    
    return(xx_df)
}