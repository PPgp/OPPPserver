#' @title Graduating Population from Five- to Single-Year Ages
#' 
#' @description The function disaggregates the input population from 5-year to 1-year ages.
#' 
#' @param pop Data table with population data aggregated to 5-year age groups. 
#'     Columns are the same as in \code{\link{get_wpp_pop}}, i.e. 
#'     \code{age}, \code{popF} and \code{popM}.
#' @param method Method to use in the \code{\link[DemoTools]{graduate}} function of the \pkg{DemoTools} package.
#' 
#' @return Data table with the same columns as in the \code{pop} input argument, but 
#'     with age groups (rows) disaggregated into single-year ages. 
#'     
#' @details The function applies the \pkg{DemoTools} function \code{\link[DemoTools]{graduate}}
#'     on both population columns of the input argument \code{pop} and returns the disaggregated dataset.
#' 
#' @export
#' 
#' @examples
#' # extract 5-year population of France in 2024 (default)
#' pop5 <- get_wpp_pop("France", n = 5)
#' 
#' # disaggregate into single-year of age
#' pop1 <- graduate_pop(pop5)
#' 
#' # plot the average  of the 5-year pop and showing it
#' # at the middle points of the age groups
#' ages5 <- seq(2, length = nrow(pop5), by = 5)
#' plot(ages5, pop5[, popF + popM]/5, type = "l", col = "blue",
#'     main = "Population of France in 2024",
#'     xlab = "age", ylab = "Population (in thousands)")
#' lines(pop1[, age], pop1[, popF + popM], col = "red")
#' legend("bottomleft", legend = c("original / 5", "graduated"),
#'     bty = "n", lty = 1, col = c("blue", "red"))
#' 
graduate_pop <- function(pop, method = "beers(ord)"){
    age5to1cat <- get_wpp("age5categories")
    age5 <- unique(age5to1cat[, list(agecat, age)])
    pop_res <- age5to1cat[, list(age = age1)]
    for(col in c("popM", "popF")){
        val <- DemoTools::graduate(pop[[col]], Age = age5$agecat, #AgeInt = rep(5, nrow(age5)),
                                   method = method, constrain = TRUE, OAG = TRUE)
        pop_res[[col]] <- val
    }
    return(pop_res)
}
