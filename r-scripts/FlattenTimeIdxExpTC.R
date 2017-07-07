#' Flatten the restricted time index given as list (per axial location) 
#'
#' Some ridiculous calculations are required when flattening the time indices
#' per axial location into a single long vector because the PCA results is 
#' missing the temperature output at 0 [s], such that the total length of time
#' points per axial location is 9'999 instead of 10'000.
#'
#' The output was excluded because original formulation of PCA was correlation
#' PCA, and correlation matrix is constructed by normalizing the elements with
#' standard deviation. At time 0 [s] the standard deviation is zero thus the 
#' data was excluded altogether.
#'
#' @param time_idx_list List of numeric vector. The list of indices indicating
#'  restricted time points per axial location
#' @param n_time Numeric. The length of time points per axial location
#' @param excl_idx Numeric. Index to exclude, but must be contiguous
#' @return
FlattenTimeIdxExpTC <- function(time_idx_list, n_time, excl_idx)
{
    time_idx_vec <- c()
    for (i in 1:length(time_idx_list))
    {
        time_idx_vec <- c(time_idx_vec, 
                          (time_idx_list[[i]][-excl_idx] - length(excl_idx)) + (i - 1) * (n_time - length(excl_idx)))
    }

    return(time_idx_vec)
}

    