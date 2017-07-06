#' Get restricted time index for Carryover output (CO <= arg [kg])
#'
#' In the series, all liquid carryover measurement saturated the tank before 
#' the end of the transient at \code{max_co} [kg], the capacity of the tank.
#' As such the calibration will be conducted with respect to part of transient
#' 
#' @param exp_data Matrix. Time-Experimental Data pairs, 1st column is the time
#' @param max_co Numeric. Maximum liquid carryover in [kg] measured in the tank
#' @return List of vectors, $idx is the time indices while $pts is the time
#'  points taken from the experimental data table with carryover measurement is 
#'  less than 10.0 [kg], inclusive one time step with measurement of 10 [kg].
GetTimeExpCO <- function(exp_data, max_co = 10.0)
{
    idx_end <- sum(exp_data[, 2] < max_co) + 1
    output <- list(
        idx = 1:idx_end,
        pts = exp_data[1:idx_end, 1]
    )

    return(output)
}
