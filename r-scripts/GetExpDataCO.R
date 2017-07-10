#' Helper function to get the vector of CO experimental data points 
#'
#' The experimental data is directly converted to [Pa] from [bar].
#'
#' @param trc_data List. List of compiled results from TRACE training runs,
#'  either PC GP or Bias GP.
#' @param max_co Numeric. The maximum carryover tank capacity 
#'  (default = 10.0 [kg])
#' @return Numeric vector. The vector of DP measurement data points in all 
#'  axial positions (bottom, middle, upper, total)
GetExpDataCO <- function(trc_data, max_co = 10.0)
{
    exp_data <- trc_data$exp_data[[3]]
    res_time <- GetTimeExpCO(exp_data, max_co = max_co)
    
    return(exp_data[res_time$idx, 2])
}
