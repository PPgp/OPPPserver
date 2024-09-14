
#' @title WPP countries
#'
#' @description Returns a vector of country names included in the WPP data.
#' @param include_aggregates Logical. If \code{TRUE} the set of countries includes
#'      aggregated regions, such as continents or other aggregations of countries.
#'      If \code{FALSE} (default) only countries are returned.
#' @param sort Logical. If \code{TRUE} results are sorted alphabetically.
#'
#' @return Character vector of names of countries and optionally aggregations.
#' @details The function uses the \code{\link[wpp2022]{UNlocations}} dataset of the
#'      \pkg{wpp2022} package. By default it selects locations where \code{location_type} is 4.
#'      If \code{include_aggregates} is \code{TRUE} it returns all locations from this dataset.
#' @export
#'
#' @examples
#' # Countries that start with "E"
#' Ecountries <- grep("^E", get_wpp_countries(), value = TRUE)
#' print(Ecountries)
#' # Europe should not be there
#' "Europe" %in% Ecountries # should be FALSE
#'
#' # include aggregations
#' Ecountries2 <- grep("^E", get_wpp_countries(include_aggregates = TRUE),
#'     value = TRUE)
#' print(Ecountries2)
#' "Europe" %in% Ecountries2 # should be TRUE as Europe is included
#'
get_wpp_countries <- function(include_aggregates = FALSE, sort = TRUE) {
    locations <- get_wpp("UNlocations")
    locations <- if(include_aggregates) locations[, "name"] else locations[locations$location_type == 4, "name"]
    if(sort) locations <- sort(locations)
    return(locations)
}

#' @title WPP regions
#'
#' @description Returns a vector of names of aggregated regions, such as continents or 
#'     other aggregations of countries, included in the WPP data.
#' @param sort Logical. If \code{TRUE} results are sorted alphabetically.
#'
#' @return Character vector of names of regions.
#' @details The function uses the \code{\link[wpp2022]{UNlocations}} dataset of the
#'      \pkg{wpp2022} package. It selects locations where \code{location_type} is other than four,
#'      as four is the code for countries.
#' @export
#'
#' @examples
#' regions <- get_wpp_regions()
#' head(regions)
#' tail(regions)
#'
get_wpp_regions <- function(sort = TRUE) {
    locations <- get_wpp("UNlocations")
    locations <- locations[locations$location_type != 4, "name"]
    if(sort) locations <- sort(locations)
    return(locations)
}


#' @title Country-specific WPP population
#'
#' @description Given a specific country and year, the function returns a sex- and
#'      age-specific dataset of population counts from the UN's World Population Prospects.
#'
#' @param country Name of country.
#' @param year Integer value specifying the year for which data should be extracted.
#'
#' @return \code{data.table} object with female and male population counts (in thousands) by age.
#'      It contains columns \code{age}, \code{popF} (female), \code{popM} (male).
#'
#' @details The data is extracted either from the dataset \code{\link[wpp2022]{popAge1dt}} or
#'      \code{\link[wpp2022]{popprojAge1dt}}, depending on the specified \code{year}.
#' @export
#'
#' @examples
#' # Retrieve and plot population of Spain for 2023 and 1960
#' spain_pop <- get_wpp_pop("Spain") # default of 2023
#' plot(spain_pop[, age], spain_pop[, popF], type = "l", col = "red",
#'     main = "Population of Spain in 2023 and 1960",
#'     xlab = "age", ylab = "Population (in thousands)")
#' lines(spain_pop[, age], spain_pop[, popM], col = "blue")
#'
#' spain_pop60 <- get_wpp_pop("Spain", year = 1960)
#' lines(spain_pop60[, age], spain_pop60[, popF], col = "red", lty = 2)
#' lines(spain_pop60[, age], spain_pop60[, popM], col = "blue", lty = 2)
#' legend("bottomleft", legend = c("female 2023", "male 2023", "female 1960", "male 1960"),
#'       bty = "n", lty = c(1, 1, 2, 2), col = rep(c("red", "blue"), 2))
#'
get_wpp_pop <- function(country, year = 2023){
    name <- NULL # to satisfy R check
    yr <- year # need to rename because a collision with the column name "year"
    pop <- get_wpp("popAge1dt") # load observed data
    if(nrow(pop[year == yr]) == 0)
        pop <- get_wpp("popprojAge1dt") # load projected data
    pop_res <- pop[name == country & as.integer(year) == yr, c("age", "popF", "popM"), with = FALSE]
    if(nrow(pop_res) == 0) stop("Either ", country, " or year ", yr, " not available in the WPP data.")
    return(pop_res)
}

#' @title Country-specific WPP total fertility rate
#'
#' @description The function returns the total fertility rate over time for a specific country
#'       from the UN's World Population Prospects.
#'
#' @param country Name of country.
#'
#' @return \code{data.table} object with columns \code{year} and \code{tfr}.
#'
#' @details The data is extracted by combining datasets \code{\link[wpp2022]{tfr1dt}} and
#'      \code{\link[wpp2022]{tfrproj1dt}}.
#' @export
#'
#' @examples
#' # Compare TFR of Niger and Brazil
#' tfr_niger <- get_wpp_tfr("Niger")
#' tfr_brazil <- get_wpp_tfr("Brazil")
#' plot(tfr_niger[, year], tfr_niger[, tfr], type = "l", ylim = c(0,8),
#'     main = "TFR of Niger and Brazil", xlab = "", ylab = "TFR")
#' lines(tfr_brazil[, year], tfr_brazil[, tfr], col = "blue")
#' legend("bottomleft", legend = c("Niger", "Brazil"), col = c("black", "blue"),
#'     lty = 1, bty = "n")
#'
get_wpp_tfr <- function(country){
    tfr_est <- get_wpp("tfr1dt") # load observed data
    tfr_proj <- get_wpp("tfrproj1dt") # load projected data
    tfr <- rbind(tfr_est[tfr_est$name == country, c("year", "tfr"), with = FALSE],
                 tfr_proj[tfr_proj$name == country, c("year", "tfr"), with = FALSE])
    if(nrow(tfr) == 0)
        stop("Country ", country, " not available in the WPP data.")
    return(tfr)
}

#' @title Country-specific WPP life expectancy at birth
#'
#' @description The function returns the life expectancy at birth for male and female over time for a specific country
#'       from the UN's World Population Prospects.
#'
#' @param country Name of country.
#'
#' @return \code{data.table} object with columns \code{year}, \code{e0M} and \code{e0F}.
#'
#' @details The data is extracted by combining datasets \code{\link[wpp2022]{e01dt}} and
#'      \code{\link[wpp2022]{e0proj1dt}}.
#' @export
#'
#' @examples
#' # Compare female e0 of Japan and Norway
#' e0_norway <- get_wpp_e0("Norway")
#' e0_japan <- get_wpp_e0("Japan")
#' plot(e0_japan[, year], e0_japan[, e0F], type = "l", ylim = c(50, 100),
#'     main = "Female e0 of Japan and Norway", xlab = "", ylab = "e0")
#' lines(e0_norway[, year], e0_norway[, e0F], col = "blue")
#' legend("bottomright", legend = c("Japan", "Norway"), col = c("black", "blue"),
#'     lty = 1, bty = "n")
#'
get_wpp_e0 <- function(country){
    e0_est <- get_wpp("e01dt") # load observed data
    e0_proj <- get_wpp("e0proj1dt") # load projected data
    e0 <- rbind(e0_est[e0_est$name == country, c("year", "e0M", "e0F"), with = FALSE],
                 e0_proj[e0_proj$name == country, c("year", "e0M", "e0F"), with = FALSE])
    if(nrow(e0) == 0)
        stop("Country ", country, " not available in the WPP data.")
    return(e0)
}

#' @title Country-specific WPP net migration counts
#'
#' @description The function returns net migration totals over time for a specific country
#'       from the UN's World Population Prospects.
#'
#' @param country Name of country.
#'
#' @return \code{data.table} object with columns \code{year} and \code{mig}.
#'
#' @details The data corresponds to the dataset \code{\link[wpp2022]{migration1dt}}.
#' @export
#'
#' @examples
#' # Compare migration of three countries
#' country.colors <- list(Germany = "black", China = "orange", 
#'                         `United States of America` = "red")
#' plot(0,1, type = "n", ylim = c(-1000, 2000), xlim = c(1950, 2100), 
#'     main = "Net Migration", ylab = "Migration/1000", xlab = "")
#' for(cntry in names(country.colors)){
#'     mig <- get_wpp_mig(cntry)
#'     lines(mig[, year], mig[, mig], col = country.colors[[cntry]])
#' }
#' abline(h = 0, col = "grey")
#' legend("topright", legend = names(country.colors), col = unlist(country.colors),
#'     lty = 1, bty = "n")
#'
get_wpp_mig <- function(country){
    mig_all <- get_wpp("migration1dt") # load observed and projected data
    mig <- mig_all[mig_all$name == country, c("year", "mig"), with = FALSE]
    if(nrow(mig) == 0)
        stop("Country ", country, " not available in the WPP data.")
    return(mig)
}



#' @title Getter function for the WPP datasets.
#' @description Since some of the datasets are big, this function caches already used datasets.
#' @keywords internal
#' @export
get_wpp <- memoise::memoise(function(dataset_name, package = "wpp2022"){
    data(list = dataset_name, package = package)
    return(get(dataset_name))
}, cache = cachem::cache_mem(max_size = 1024 * 1024^2,
                             max_age = 3600*3)) # 1G for 3 hours

get_country_code <- function(country){
    # return the UN code for given country
    locations <- get_wpp("UNlocations")
    column <- if(is.numeric(country)) "country_code" else "name"
    return(locations[locations[[column]] == country, "country_code"])
}

get_wpp_indicator_multiple_years <- function(indicator_est, indicator_proj = NULL, un_code, end_year = NULL, ...){
    country_code <- NULL # to satisfy CRAN check
    all_wpp <- get_wpp(indicator_est, ...)[country_code == un_code] # load observed data
    if(!is.null(indicator_proj) && (is.null(end_year) || end_year > all_wpp[, max(year)])) # add projected data if needed
        all_wpp <- rbind(all_wpp, get_wpp(indicator_proj, ...)[country_code == un_code],
                             fill = TRUE)[, year := as.integer(year)]
    if(!is.null(end_year)) all_wpp <- all_wpp[year <= end_year] # end_year should be the latest year in the data
    return(all_wpp)
}
get_wpp_pop_by_age_multiple_years <- function(un_code, end_year = NULL){
    return(get_wpp_indicator_multiple_years("popAge1dt", "popprojAge1dt", un_code = un_code, end_year = end_year))
}

get_wpp_pop_multiple_years <- function(un_code, end_year = NULL){
    return(get_wpp_indicator_multiple_years("pop1dt", "popproj1dt", un_code = un_code, end_year = end_year))
}
