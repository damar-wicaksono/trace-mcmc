#' Helper function to get the vector of DP experimental data points 
#'
#' The experimental data is directly converted to [Pa] from [bar].
#'
#' @param trc_data List. List of compiled results from TRACE training runs,
#'  either PC GP or Bias GP.
#' @return Numeric vector. The vector of DP measurement data points in all 
#'  axial positions (bottom, middle, upper, total)
GetExpDataDP <- function(trc_data)
{
    exp_data <- trc_data$exp_data[[2]]
    exp_vec <- c()
    for (j in 1:4)
    {
        exp_vec <- c(exp_vec, exp_data[, j+1])
    }
    exp_vec <- exp_vec * 1e5    # Convert to [Pa]
}