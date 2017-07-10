#' Helper function to get the vector of indices to subset the full CO output
#'
#' The available measurement data is maxed at 10 [kg].
#' As such compared the TRACE prediction at time points before.
#'
#' @param trc_data List. List of compiled results from TRACE model training 
#'  either PC GP or Bias GP runs
#' @param d_t Numeric. The time step size in the uniform time grid of full 
#'  output space (default = 0.1 [s])
#' @param max_co Numeric. The maximum carryover tank capacity 
#'  (default = 10.0 [kg])
#' @return Numeric vector. The indices to subset the full model output of CO
#'  up to 10.0 [kg] liquid carryover
GetTimeIdxCO <- function(trc_data, d_t = 0.1, max_co = 10.0)
{
    exp_data <- trc_data_bias$exp_data[[3]]
    
    res_time <- GetTimeExpCO(exp_data, max_co = max_co)
    time_idx <- GetTimeIdxTrc(res_time$pts, d_t = d_t)
    
    return(time_idx)
}
