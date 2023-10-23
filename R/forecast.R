
#' @title Run Population Forecast
#' @description Runs forecast and returns various datasets with results.
#'
#' @param country Name of country.
#' @param start_year Start of the projection. It is the last observed year; the first projected year
#'     is \code{start_year + 1}.
#' @param end_year The last year of the projection.
#' @param pop Data table with the initial population data at \code{start_year}.
#'     Columns are the same as in \code{\link{get_wpp_pop}}. By default, the WPP data is used.
#' @param tfr Data table containing the TFR scenario. It has the same structure as
#'     the dataset returned by \code{\link{get_wpp_tfr}}. By default the UN data is used.
#' @param units Scale of the output data. By default a simulation results are given in thousands.
#'     If for example one would like to see results in millions, set \code{unit = 1e6}.
#'
#' @return List with elements:
#' \describe{
#'     \item{prediction}{\code{\link[bayesPop]{bayesPop.prediction}} object.
#'          (This is for debugging and will be removed in the future.)}
#'     \item{population_by_age_and_sex}{Data table with age- and sex-specific population by time
#'          (columns \code{year}, \code{age}, \code{popM}, \code{popF}).}
#'      \item{population_by_broad_age_group}{Data table where population is grouped by broader age groups.
#'          (columns \code{year}, \code{age}, \code{pop}, \code{pop_percent}).}
#'      \item{population_by_time}{Data table where with total population as well as by broader age groups.
#'          It also contains the UN median and the 95\% probability intervals 
#'          (columns \code{year}, \code{age}, \code{pop}, \code{un_pop_median}, 
#'          \code{un_pop_95low}, \code{un_pop_95high}).}
#'      \item{tfr_by_time}{Data table with the total fertility forecast and the UN median and the 
#'          95\% probability intervals (columns \code{year}, \code{tfr}, \code{un_tfr_median}, 
#'          \code{un_tfr_95low}, \code{un_tfr_95high}).}
#'      \item{annual_growth_rate}{Data table with the forecast of annual growth rate as a percentage,
#'              computed as \eqn{log(P_t/P_{t-1})*100}, for the same age groups as in \code{population_by_time}
#'              (columns \code{year}, \code{age}, \code{growth_rate}).}
#'      \item{birth_count_rates}{Data table with the forecast of total births (column \code{births}),
#'              crude birth rate (\code{cbr}) and the UN median and 95\% probability intervals 
#'              of the respective quantities
#'              (columns \code{un_births_median}, \code{un_births_95low}, \code{un_births_95high}, 
#'              \code{un_cbr_median}, \code{un_cbr_95low}, \code{un_cbr_95high})}
#' }
#' @details For now the function only runs when \code{start_year} is 2021 or smaller.
#'     Note that for now each run of this function creates a new temporary directory.
#'     It can be removed by running \code{remove_forecast()} as in the example below.
#'
#' @examples
#' # Change the TFR of Niger to stay at the 4.0 level from 2054 onwards
#' my_tfr <- get_wpp_tfr("Niger")
#' my_tfr[year >= 2054, tfr := 4]
#'
#' forecast <- run_forecast("Niger", start_year = 2021, tfr = my_tfr)
#'
#' # plot population and fertility
#' ################################
#' par(mfrow = c(1,2))
#' 
#' # extract total population
#' respop <- forecast$population_by_time[age == "0+"]
#'
#' # the scenario population should be larger than UN starting in 1954
#' plot(respop$year, respop$pop, type = "l", col = "red",
#'     main = "Niger population", ylab = "Population", xlab = "")
#' lines(respop$year, respop$un_pop_median, col = "blue")
#' lines(respop$year, respop$un_pop_95low, col = "blue", lty = 2)
#' lines(respop$year, respop$un_pop_95high, col = "blue", lty = 2)
#' legend("topleft", legend = c("my scenario", "UN", "UN 95% PI"), col = c("red", "blue", "blue"),
#'     lty = c(1, 1, 2), bty = "n")
#'     
#' # extract fertility results
#' restfr <- forecast$tfr_by_time
#' plot(restfr$year, restfr$tfr, type = "l", col = "red", ylim = c(0, 8),
#'     main = "Niger TFR", ylab = "Total fertility rate", xlab = "")
#' lines(restfr$year, restfr$un_tfr_median, col = "blue")
#' lines(restfr$year, restfr$un_tfr_95low, col = "blue", lty = 2)
#' lines(restfr$year, restfr$un_tfr_95high, col = "blue", lty = 2)
#' legend("bottomleft", legend = c("my scenario", "UN", "UN 95% PI"), col = c("red", "blue", "blue"),
#'     lty = c(1, 1, 2), bty = "n")
#' 
#'
#' # remove prediction directory
#' remove_forecast(forecast)
#'
#' @export
#'

run_forecast <- function(country, start_year = 2021, end_year = 2100,
                         pop = NULL, tfr = NULL, units = 1000){
    # For now it only works if start_year <= 2021!

    code <- get_country_code(country)

    # create input files for pop.predict
    pop_files <- prepare_pop(pop, code, start_year)
    tfr_file <- prepare_tfr(tfr, country, code, start_year)

    # simulation directory
    sim_dir <- tempdir()

    # launch simulation
    ###################
    pred <- pop.predict(end.year = end_year, present.year = start_year,
                        wpp.year = 2022, output.dir = sim_dir,
                        countries = code, annual = TRUE,
                        inputs = list(
                            popM = pop_files["male"], popF = pop_files["female"],
                            tfr.file = tfr_file
                        ),
                        nr.traj = 1, keep.vital.events = TRUE,
                        fixed.mx = TRUE,
                        replace.output = TRUE
                        )

    # In order not to copy the pred object every time we need it in a function,
    # we put it into an environment that acts as a pointer
    pred_env <- new.env()
    pred_env[["prediction"]] <- pred
    
    # collect results
    #################
    pop_by_age_sex <- extract_pop_by_age_sex(pred_env, units = units)
    pop_by_broad_age <- get_pop_by_broad_age(pop_by_age_sex)
    pop_by_time <- get_pop_by_time(pop_by_age_sex, code)
    tfr_by_time <- extract_tfr_by_time(pred_env)
    growth_rate <- get_annual_growth_rate(pop_by_time)
    births <- extract_births(pred_env, units = units)

    # cleanup
    ##########
    unlink(pop_files)

    # For now include the whole prediction object in the return value (easier for testing).
    # It should be cleaned up via the remove_forecast() function.
    # Be aware that each call to run_forecast creates a new temp directory!
    # Later we should uncomment this line and not return the "pred" object.
    #unlink(sim_dir, recursive = TRUE)

    return(list(prediction = pred,
                population_by_age_and_sex = pop_by_age_sex,
                population_by_broad_age_group = pop_by_broad_age,
                population_by_time = pop_by_time,
                tfr_by_time = tfr_by_time,
                annual_growth_rate = growth_rate,
                birth_count_rates = births
                )
           )
}


prepare_pop <- function(pop, un_code, start_year){
    # Stores population into two files (one for each sex)
    # that are compatible with bayesPop
    # (items popM & popF of "inputs" argument in ?pop.predict)
    # Here we will replace the WPP data of the particular year
    # with the new data.
    if(is.null(pop)) return(c(male = NULL, female = NULL))

    country_code <- i.popM <- i.popF <- NULL # to satisfy CRAN check

    # load UN data
    all_wpp_pop <- get_wpp_pop_by_age_multiple_years(un_code, start_year)

    # replace data at the start_year with the user-provided data and add country_code
    popext <- cbind(pop, year = start_year)
    all_wpp_pop[popext, `:=`(popM = i.popM, popF = i.popF),
                on = c("year", "age")][, country_code := un_code]

    # convert to wide format
    male_pop <- dcast(all_wpp_pop, country_code + age ~ year, value.var = "popM")
    female_pop <- dcast(all_wpp_pop, country_code + age ~ year, value.var = "popF")

    # save to temp files
    male_file <- tempfile("popM", fileext = ".txt")
    female_file <- tempfile("popF", fileext = ".txt")
    fwrite(male_pop, file = male_file, sep = "\t")
    fwrite(female_pop, file = female_file, sep = "\t")
    return(c(male = male_file, female = female_file))
}

prepare_tfr <- function(tfr, country, un_code, start_year){
    # Stores TFR into a file that are compatible with bayesPop
    # (item tfr.file of "inputs" argument in ?pop.predict)
    # The WPP data for all years provided in the tfr dataset
    # are replaced by the new data.
    if(is.null(tfr)) return(NULL)

    country_code <- Trajectory <- i.tfr <- NULL # to satisfy CRAN check

    all_wpp_tfr <- get_wpp_tfr(country)[, country_code := un_code]

    # replace data with the user-provided data
    all_wpp_tfr[tfr, tfr := i.tfr, on = "year"]

    # as this represents a TFR projection, start_year should be the smallest year in the dataset
    #all_wpp_tfr <- all_wpp_tfr[year >= start_year]

    # prepare for storing into a csv trajectory file
    setnames(all_wpp_tfr[, Trajectory := 1],
             c("country_code", "year", "tfr"),
             c("LocID", "Year", "TF"))

    # save to a temp file
    tfr_file <- tempfile("tfr", fileext = ".csv")
    fwrite(all_wpp_tfr, file = tfr_file, sep = ",")
    return(tfr_file)
}

extract_pop_by_age_sex <- function(env, units = 1000){
    age <- age_to_100 <- pop <- NULL # to satisfy CRAN check
    pop_dt_sx <- list()
    country <- env$prediction$countries$code
    for(sx in c("M", "F")){
        # extract observed pop
        obs_df <- get.pop.exba(paste0("P", country, "_", sx, "{}"),
                               env$prediction, observed = TRUE)
        if(units != 1000)
            obs_df <- obs_df*1000 / units

        # make it data.table and convert to long format
        obs_dt <- cbind(data.table(age = 0:(nrow(obs_df)-1)), data.table(obs_df))
        obs_dt_long <- melt(obs_dt, id.vars = "age", variable.name = "year",
                value.name = paste0("pop", sx), variable.factor = FALSE)

        # extract the projected median matrix
        pred_df <- env$prediction[[paste0("quantiles", sx, "age")]][1, , "0.5",]
        if(units != 1000)
            pred_df <- pred_df*1000 / units

        # make it data.table
        pred_dt <- cbind(data.table(age = rownames(pred_df)), data.table(pred_df))

        # convert to long format
        pred_dt_long <- melt(pred_dt, id.vars = "age", variable.name = "year",
                                value.name = "pop", variable.factor = FALSE)

        # remove years that are already in the observed data (there is an overlap of the start_year)
        pred_dt_long <- pred_dt_long[! year %in% obs_dt_long[, year]]

        # since projections are up to age 130, aggregate high ages into 100+
        pred_dt_long[, age := as.integer(age)]
        pred_dt_long[, age_to_100 := ifelse(age <= 100, age, 100)]
        pred_dt_to100 <- pred_dt_long[, list(pop = sum(pop)),
                                      by = c("year", "age_to_100")]

        # rename the indicator and age columns
        setnames(pred_dt_to100, c("pop", "age_to_100"),
                                c(paste0("pop", sx), "age"))

        # join observed and predicted
        pop_dt_sx[[sx]] <- rbind(obs_dt_long, pred_dt_to100)[, year := as.integer(year)]
    }
    # merge both sexes together
    return(merge(pop_dt_sx$M, pop_dt_sx$F, by = c("year", "age")))
}

get_pop_by_broad_age <- function(pop_by_age_sex,
                                 age_groups = c(0, 20, 40, 60, 150),
                                 compute_percent = TRUE){
    age <- pop <- year <- popM <- popF <- age_group <- sum_pop <- pop_percent <- NULL # to satisfy CRAN check

    # sum population over sexes
    dt <- pop_by_age_sex[, list(year, age, pop = popM + popF)]

    # identify broad age groups
    dt[, age_group := cut(age, age_groups, labels = FALSE,
                          include.lowest = TRUE, right = FALSE)]
    # sum over broad age groups
    dt_broad <- dt[, list(pop = sum(pop)), by = c("year", "age_group")]

    # compute percent
    if(compute_percent)
        dt_broad[, sum_pop := sum(pop), by = "year"][
            , pop_percent := pop/sum_pop * 100][, sum_pop := NULL]

    # compose labels for the broad age groups
    dt_broad[age_group < (length(age_groups) - 1), # for all closed age groups
             age := paste(age_groups[age_group],
                          age_groups[age_group + 1] - 1, sep = "-")]
    dt_broad[age_group == length(age_groups) - 1,  # open age group
             age := paste0(age_groups[age_group], "+")
             ][, age_group := NULL]

    # re-order columns
    setcolorder(dt_broad, c("year", "age"))
    return(dt_broad)
}

get_pop_by_time <- function(pop_by_age_sex, un_code) {
    age <- un_pop_median <- i.pop <- i.pop_95l <- i.pop_95u <- NULL # to satisfy CRAN check
    # totals
    dt <- get_pop_by_broad_age(pop_by_age_sex, age_groups = c(0, 150),
                               compute_percent = FALSE)
    # standard broad age groups
    dt <- rbind(dt, get_pop_by_broad_age(pop_by_age_sex,
                                         compute_percent = FALSE))
    # 65+
    dt <- rbind(dt, get_pop_by_broad_age(pop_by_age_sex[age >= 65],
                                         age_groups = c(65, 150),
                                         compute_percent = FALSE))

    # Add UN prediction median
    unpop <- get_wpp_pop_by_age_multiple_years(un_code)
    unpop_groups <- get_pop_by_broad_age(unpop, age_groups = c(0, 150),
                                         compute_percent = FALSE)
    unpop_groups <- rbind(unpop_groups,
                          get_pop_by_broad_age(unpop, compute_percent = FALSE))
    unpop_groups <- rbind(unpop_groups,
                          get_pop_by_broad_age(unpop[age >= 65],
                                               age_groups = c(65, 150),
                                               compute_percent = FALSE))
    dt[unpop_groups, un_pop_median := i.pop, on = c("year", "age")]

    # add prediction intervals for total pop
    untotpop <- get_wpp_pop_multiple_years(un_code)[, age := "0+"]
    # merge with dt for the 0+ age group
    dt[untotpop, `:=`(un_pop_95low = i.pop_95l, un_pop_95high = i.pop_95u), 
       on = c("year", "age")]
    
    # TODO: add intervals for broader age groups 
    # (need to find them on the UN site and put somewhere)
    return(dt)
}

extract_tfr_by_time <- function(env) {
    tfr <- tfr_95l <- tfr_95u <- un_tfr_median <- NULL
    # extract observed and predicted data
    expression <- paste0("F", env$prediction$countries$code)
    tfr_data <- c(get.pop.ex(expression, env$prediction, observed = TRUE),
            get.pop.ex(expression, env$prediction, observed = FALSE)[-1])
    
    # convert to long format
    tfr_long <- data.table(year = as.integer(names(tfr_data)), tfr = tfr_data)
    
    # extract UN median and PI intervals
    untfr <- get_wpp_indicator_multiple_years("tfr1dt", "tfrproj1dt", 
                                              un_code = env$prediction$countries$code)
    
    # merge together
    dt <- merge(tfr_long, untfr[, list(year, un_tfr_median = tfr, 
                                    un_tfr_95low = tfr_95l, un_tfr_95high = tfr_95u)],
                all = TRUE, by = "year")
    return(dt[!is.na(un_tfr_median)])
}

get_annual_growth_rate <- function(pop){
    age <- year_lag <- pop_lag <- i.pop <- NULL
    popdt <- copy(pop)[, year_lag := year - 1]
    popdt[popdt, pop_lag := i.pop, on = c(year_lag = "year", age = "age")]
    return(popdt[!is.na(pop_lag), list(year, age, growth_rate = log(pop/pop_lag)*100)])
}

extract_births <- function(env, units = 1000){
    births <- births_95l <- births_95u <- cbr <- cbr_95l <- cbr_95u <- NULL
    # extract observed and predicted counts
    uncode <- env$prediction$countries$code
    expression <- paste0("B", uncode)
    birth_count <- c(get.pop.ex(expression, env$prediction, observed = TRUE),
                  get.pop.ex(expression, env$prediction, observed = FALSE)[-1])
    if(units != 1000)
        birth_count <- birth_count*1000 / units
    
    # extract observed and predicted rates
    expression <- paste0("1000 * B", uncode, "/mid.period(P", uncode, ")")
    birth_rate <- c(get.pop.ex(expression, env$prediction, observed = TRUE),
                     get.pop.ex(expression, env$prediction, observed = FALSE)[-1])
    
    # convert to long format
    births_long <- data.table(year = as.integer(names(birth_count)), 
                              births = birth_count, cbr = birth_rate)
    
    # extract UN values
    unbirths <- get_wpp_indicator_multiple_years("births1dt", "birthsproj1dt", package = "wpp2022extra",
                                              un_code = uncode)
    uncbr <- get_wpp_indicator_multiple_years("cbr1dt", "cbrproj1dt", package = "wpp2022extra",
                                              un_code = uncode)
    # merge together
    dt <- merge(
            merge(births_long, 
                unbirths[, list(year, un_births_median = births, 
                          un_births_95low = births_95l, un_births_95high = births_95u)],
                all = TRUE, by = "year"
                ),
            uncbr[, list(year, un_cbr_median = cbr, 
                        un_cbr_95low = cbr_95l, un_cbr_95high = cbr_95u)],
                all = TRUE, by = "year"
            )
    return(dt[!is.na(births)])
}


#' @export
remove_forecast <- function(forecast)
    unlink(forecast$base.directory, recursive = TRUE)
