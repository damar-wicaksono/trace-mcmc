#' Get restricted time index for Carryover output (CO <= arg [kg])
#'
#' \code{GetTimeIdxCO} takes the last time index in a table of time-data pairs
# of liquid carryover measurement in FEBA Test Series I.
# In the series, all liquid carryover measurement saturated the tank before 
# the end of the transient at \code{max_co} [kg], the capacity of the tank.
# As such the calibration will be conducted with respect to part of transiet
# 
# @param exp_data Matrix. Time-Experimental Data pairs, 1st column is the time
# @param max_co Numeric. Maximum liquid carryover in [kg] measured in the tank
# @return Vector of time indices in the experimental data table whose carryover
#   measurement is less than 10.0 [kg], and inclusive one time step with 
#   measurement of 10 [kg].
GetTimeIdxExpCO <- function(exp_data, max_co = 10.0)
{
    idx_end <- sum(exp_data[,2] < max_co) + 1
    return(1:idx_end)
}
