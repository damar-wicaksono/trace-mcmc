#' Create a single element dataframe with consistent naming for PC Scores GP 
#'
#' The input to GP metamodel of PC scores has to be a dataframe specifically 
#' having the same number of columns (variables) AND the same name. 
#' \code{CreateInputPCScores} is a helper function to construct a single
#' element dataframe with names as parameter, with the input values coming 
#' from the MCMC sampler.
#' 
#' @param xx Numeric vector. Input values for GP metamodel of PC scores
#' @param str_names Character vector. Names of the input parameters
#' @return Dataframe. A single element dataframe to be used in GP metamodel of
#'  PC scores
CreateInputPCScores <- function(xx, str_names)
{
    xx_df <- data.frame(matrix(0, ncol = length(xx), nrow = 0))
    colnames(xx_df) <- str_names
    xx_df[1,] <- xx
    return(xx_df)
}