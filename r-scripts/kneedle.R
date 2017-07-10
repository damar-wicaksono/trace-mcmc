#' Calculate the knee point from a curve
#'
#' Knee point is an important feature of the reflood curves dataset.
#' \code{kneedle} is an implementation of the "kneedle" algorithm
#' Reference:
#'   V. Satopaa et al., "Finding a "Kneedle" in a Haystack: Detecting Knee 
#'   Points in System Behavior," In. Proc. the 31st International Conference 
#'   on Distributed Computing Systems Workshop (ICDCSW 11), Washington, 
#'   pp. 166 - 177, 2011.
#' Use findpeaks.R
#' 
#' @param x A vector of float. The argument.
#' @param y A vector of float. The function values
#' @return the argument value the gives the knee point
kneedle <- function(x, y) {
    # An implementation of the "kneedle" algorithm.
    
    # The first data always create problem, remove them from the global max.
    # selection. Here we only consider part of the transient after the max.
    x = x[which.max(y):length(x)]
    y = y[which.max(y):length(y)]
    
    # Normalized the x and y data
    normalized_x <- (x - min(x)) / (max(x) - min(x))
    normalized_y <- (y - min(y)) / (max(y) - min(y))
    
    # Invert the data to have a knee curve
    normalized_y_inverted <- -1.0 * normalized_y + 1.0
    
    # Calculate the difference between the curve and the y = x curve
    difference_data <- normalized_y_inverted - normalized_x
    
    # Calculate the global maximum of the difference_data
    # This defines the knee point as a landmark
    
    # Vector of local maxima indices
    local_max_index <- findpeaks(difference_data)[[1]]
    # Global maximum for difference_data
    global_max_index <- which.max(difference_data[local_max_index])
    #global_max_index <- sort(difference_data[local_max_index], index.return=T, 
    #                         decreasing=T)$ix[2]
    global_max_x_normalized <- normalized_x[local_max_index[global_max_index]]
    
    # Renormalized x
    global_max_x <- min(x) + global_max_x_normalized * ( max(x) - min(x) )
    
    return(global_max_x)
}
