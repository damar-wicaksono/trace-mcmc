#' Half Cauchy Density (x < 0 is forced to be 0)
#'
#' The density of half Cauchy is twice the density of standard Cauchy because
#' the area under the curve has to be equal to 1.0
#' The default scale parameter value is re-defined to be 25.0
#' 
#' @param x Numeric vector. The vector of quantile
#' @param location Numeric. The location parameter
#' @param scale  Numeric. The scale parameter
#' @param log Boolean. Log scale computation
#' @return the density of the half Cauchy distribution
dhalfcauchy <- function(x, location = 0, scale = 25, log = FALSE)
{
    if (all(x < 0))
    {
        if (log) {return(log(0))} else {return(0)}
    } else
    {   
        x <- x[x >= 0]  # Exclude negative value
        
        return(2 * dcauchy(x, location = location, scale = scale, log = log))
    }
}