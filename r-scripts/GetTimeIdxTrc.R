#' Get the time index for TRACE based on the experimental data time points
#' 
#' Note that the division between any of the time points and d_t has to yield
#' positive integer value (strictly speaking, whole number) 
#'
#' @param exp_time Numeric vector. Time points of the experimental data
#' @param d_t Numeric. Time step size of the TRACE simulation.
#'   It is assumed that the step size is uniform, that the time grid has been 
#'   homogenized.
#' @return Vector of time indices in the TRACE output
GetTimeIdxTrc <- function(exp_time, d_t)
{
    if (sum(exp_time %% d_t) != 0)
    {
        stop("One or more time index is not a whole number")   
    } else 
    {
        return(exp_time / d_t + 1)
    }
}
