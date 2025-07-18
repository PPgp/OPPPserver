#' @title Graduating Population from Five- to Single-Year Ages
#' 
#' @description The function disaggregates the input population from 5-year to 1-year ages.
#' 
#' @param pop Data table with population data aggregated to 5-year age groups. 
#'     The only required column is \dQuote{age} with starting values of the age groups,
#'     i.e. 0, 5, 10, \dots.  Any additional columns are considered as population 
#'     values to be disaggregated.
#' @param method Method to use in the \code{\link[DemoTools]{graduate}} function 
#'     of the \pkg{DemoTools} package.
#' @param pop_columns Character vector of columns of the \code{pop} dataset to graduate. 
#'     By default, all columns that are not called \dQuote{age} are used.
#' 
#' @return Data table with the same columns as in the \code{pop} dataset 
#'     (or defined by \code{pop_columns}), 
#'     with age groups (rows) disaggregated into single-year ages. 
#'     
#' @details The function applies the \pkg{DemoTools} function \code{\link[DemoTools]{graduate}}
#'     on all population columns of the \code{pop} dataset and returns the disaggregated dataset.
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
#' # plot the average of the 5-year pop at the middle of each age group
#' plot(pop5[, age + 2], pop5[, popF + popM]/5, type = "l", col = "blue",
#'     main = "Population of France in 2024",
#'     xlab = "age", ylab = "Population (in thousands)")
#' lines(pop1[, age], pop1[, popF + popM], col = "red")
#' legend("bottomleft", legend = c("original / 5", "graduated"),
#'     bty = "n", lty = 1, col = c("blue", "red"))
#' 
graduate_pop <- function(pop, method = "beers(ord)", pop_columns = NULL){
    age <- agecat <- age1 <- NULL
    age5to1cat <- get_wpp("age5categories")
    age5 <- unique(age5to1cat[, list(agecat, age)])
    pop_res <- age5to1cat[age1 <= max(pop[["age"]]), list(age = age1)]
    if(is.null(pop_columns)) pop_columns <- setdiff(colnames(pop), "age")
    for(col in pop_columns){
        val <- DemoTools::graduate(pop[[col]], Age = pop[["age"]], 
                                   method = method, constrain = TRUE, OAG = TRUE)
        pop_res[[col]] <- as.vector(val)
    }
    return(pop_res)
}

#' @title Collapsing Old Ages
#' @description Collapse old ages into a new open age group (OAG).
#' 
#' @param pop Data table with columns \dQuote{age} and other columns containing 
#'     age-specific population data. It can be either for single year or five year age groups.
#' @param oag_new Integer defining the new OAG.It should be smaller
#'      than the current OAG, i.e. the maximum age available in the \code{pop} dataset.
#'      If \code{pop} is defined for five year age groups,
#'     \code{oag_new} must be a multiple of 5.
#' @param pop_columns Character vector of columns of the \code{pop} dataset to collapse. 
#'     By default, all columns that are not called \dQuote{age} are used.
#'     
#' @return Data table with the same columns as in the \code{pop} dataset 
#'     (or defined by \code{pop_columns}), where all age groups equal to and older than 
#'     \code{oag_new} are aggregated.
 
#' @export
#' 
#' @examples
#' # extract population of Canada in 1970
#' pop_full <- get_wpp_pop("Canada", year = 1970)
#' 
#' # sum population by sex
#' pop_full[, pop_tot := popM + popF]
#' 
#' # collapse into OAG of 95
#' pop_reduced <- reduce_oag(pop_full, oag_new = 95)
#' 
#' tail(pop_full)
#' tail(pop_reduced)
#' 
#' # works for 5-year age groups as well
#' pop5_full <- get_wpp_pop("Canada", n = 5)
#' pop5_reduced <- reduce_oag(pop5_full, oag_new = 80)
#' 
#' tail(pop5_full)
#' tail(pop5_reduced)
#' 
reduce_oag <- function(pop, oag_new, pop_columns = NULL) {
    age <- NULL
    if(is.null(pop_columns)) pop_columns <- setdiff(colnames(pop), "age")
    if(! oag_new %in% pop[["age"]])
        stop("Wrong argument oag_new. It must be found in the age column of pop.")
    return(rbind(pop[age < oag_new], 
                     cbind(age = oag_new, 
                           pop[age >= oag_new, lapply(.SD, sum), .SDcols = pop_columns]))
            )
}

#' @title Extrapolating Open Age group
#' @description Extrapolate the open age group (OAG) into multiple ages using a standard population.
#' 
#' @param pop Data table of population counts with columns 
#'     \dQuote{age}, \dQuote{popM} and \dQuote{popF}. Ages can be either 
#'     single year or five year age groups. 
#' @param country Name of the country to use for the population standard.
#' @param year Which year should be used for extracting the standard population.
#' @param oag_new Integer defining the new OAG. It must be available in WPP population datasets,
#'     i.e. it must be smaller equal 100. If \code{pop} is defined for five year age groups,
#'     \code{oag_new} must be a multiple of 5.
#' @param n Size of age groups. It can be either 1 (default) or 5. It must correspond 
#'     to the age grouping in the \code{pop} dataset.
#'     
#' @return Data table with the same columns as in the \code{pop} dataset, 
#'     where the original last age group is dissaggregated up to the age 
#'     given by \code{oag_new}.
#'     
#' @export
#' 
#' @examples
#' # extract population of Germany in 2023 
#' pop_full <- get_wpp_pop("Germany", year = 2023)
#' 
#' # collapse the OAG into 85
#' pop_oag85 <- reduce_oag(pop_full, 85)
#' 
#' # extrapolate to 95
#' pop_oag95 <- extend_oag(pop_oag85, "Germany", year = 2023, oag_new = 95)
#' 
#' # plot results
#' plot(pop_oag85[, age], pop_oag85[, popM], type = "l", xlim = c(50, 100),
#'     ylim = range(pop_oag85[, popM], pop_full[, popM]),
#'     main = "Male Population of Germany in 2023", col = "blue",
#'     xlab = "age", ylab = "Population (in thousands)")
#' lines(pop_full[, age], pop_full[, popM], col = "black")
#' lines(pop_oag95[, age], pop_oag95[, popM], col = "red")
#' legend("bottomleft", legend = c("full", "OAG = 85", "OAG = 95"),
#'     bty = "n", lty = 1, col = c("black", "blue", "red"))
#'
#' # check the sums for pop of single year ages
#' all.equal(pop_full[, sum(popF + popM)], pop_oag85[, sum(popF + popM)])
#' all.equal(pop_full[, sum(popF + popM)], pop_oag95[, sum(popF + popM)])
#' 
#' # 5-year age groups
#' pop5_full <- get_wpp_pop("Germany", year = 2023, n = 5)
#' pop5_oag80 <- reduce_oag(pop5_full, 80)
#' pop5_oag100 <- extend_oag(pop5_oag80, "Germany", year = 2023, oag_new = 100, n = 5)
#' # check the sums
#' all.equal(pop5_full[, sum(popF + popM)], pop5_oag80[, sum(popF + popM)])
#' all.equal(pop5_full[, sum(popF + popM)], pop5_oag100[, sum(popF + popM)])
#' 

extend_oag <- function(pop, country, year = 2023, oag_new = 100, n = 1){
    age <- NULL
    oag_now <- max(pop[["age"]])
    pop_full <- get_wpp_pop(country, year, n = n)
    if(! oag_new %in% pop_full[["age"]])
        stop("Wrong argument oag_new. It must be an age available in the WPP population, i.e. <= 100.")
    if(length(setdiff(pop_full[age <= oag_now][["age"]], pop[["age"]])) > 0)
        stop(paste("Ages of the population standard do not correspond to ages in the given pop dataset. Check if the argument 'n' is correct.",
                   "\n\tStandard: ", paste(pop_full[age <= oag_now][["age"]], collapse = ","), 
                   "\n\tGiven: ", paste(pop[["age"]], collapse = ",")))
    if(oag_new == oag_now) return(pop)
    if(oag_new < oag_now) return(reduce_oag(pop, oag_new))
    pop_res <- pop_full[age <= oag_new, list(age)]
    for(col in intersect(c("popM", "popF"), colnames(pop))){
        pop_res[[col]] <- DemoTools::OPAG_simple(pop[[col]], Age = pop[["age"]],
                                                 OAnow = oag_now, StPop = pop_full[[col]],
                                                 StAge = pop_full[["age"]], OAnew = oag_new)
    }
    return(pop_res)
}