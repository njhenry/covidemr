
#' Create age group categories
#'
#' Given a set of age cutoffs starting at 0, this function deterministically
#' creates a data.table of age groupings with IDs and lower/upper bounds for
#' each age group
#'
#' @param age_cutoffs Vector of lower bounds for each age group, starting at 0
#'
#' @return Data.table containing age group categories
#'
create_age_groups <- function(age_cutoffs = seq(0, 100, by=20)){
  if(length(age_cutoffs)==0) stop("Need at least one age group!")
  age_groups <- data.table(age_group_code = 1:length(age_cutoffs))
  # Set start and end years for each age group
  age_groups[, age_group_years_start := age_cutoffs ]
  if(length(age_cutoffs)==1){
    age_groups[, age_group_years_end := as.numeric(NA)]
  } else {
    age_groups[, age_group_years_end := c(age_cutoffs[2:length(age_cutoffs)]-1, NA)]
  }
  # Add a name for each age group
  age_groups[, age_group_name := paste(age_group_years_start, 'to', age_group_years_end) ]
  age_groups[.N, age_group_name := paste(age_group_years_start, 'and above')]
  return(age_groups)
}


#' Prepare death data for Italy
#'
#' Given raw data all-cause mortality by age, sex, date, and province in Italy,
#' prepare the data using standard identifiers and grouped by year and week.
#' This is a convenience function designed specifically to work with data
#' downloaded from Istat. For more information, see the repository README.
#'
#' @param deaths_raw Data.table downloaded from IStat
#' @param age_cutoffs Vector of the starting years for each age group bin. For
#'   more information, see \link{\code{create_age_groups}}
#' @param model_years Vector of years to be included in the mortality model, up
#'   to and including the current year
#' @param first_covid_death_date R Date object designating the date of the first
#'   confirmed COVID death in Italy.
#'
#' @return Data.table of formatted death data
#'
ita_prepare_deaths <- function(
  deaths_raw, age_cutoffs, model_years, first_covid_death_date
){
  # Keep only required columns and reshape long
  id_cols <- c('PROV','CL_ETA','GE')
  val_cols <- c(outer(c('M_', 'F_'), model_years - 2000, FUN='paste0'))
  deaths_long <- data.table::melt(
    deaths_raw,
    id.vars = id_cols,
    measure.vars = val_cols,
    variable.name = 'sex_yr',
    value.name = 'deaths'
  )
  setnames(deaths_long, id_cols, c('location_code','ag_orig','day_id'))
  deaths_long[, sex_yr := as.character(sex_yr) ]
  deaths_long[, day_id := sprintf("%04d", day_id) ]

  # Drop deaths from non-observed days
  deaths_long <- deaths_long[!is.na(deaths)]

  # Add age bins
  deaths_long[, age_years := (ag_orig - 1) * 5]
  deaths_long[ ag_orig %in% c(0,1), age_years := ag_orig]
  # Add age groups
  age_groups <- create_age_groups(age_cutoffs=age_cutoffs)
  deaths_long$age_group_code <- age_groups$age_group_code[
    findInterval(deaths_long$age_years, age_groups$age_group_years_start)
  ]
  # Add sexes and years
  deaths_long[, sex := 'male']
  deaths_long[ startsWith(sex_yr, 'F'), sex := 'female']
  deaths_long[, year := as.integer(substr(sex_yr, 3, 4)) + 2000 ]
  # Add month and day, then convert to week ID
  deaths_long[, month := as.integer(substr(day_id, 1, 2))]
  deaths_long[, day := as.integer(substr(day_id, 3, 4)) ]
  deaths_long[, date_full := as.Date(sprintf('%04d-%02d-%02d', year, month, day)) ]
  # Weeks are defined as starting on January 1st (Jan 8 starts week 2, etc.)
  deaths_long[, week := as.integer(ceiling(
    as.integer(strftime(date_full, format='%j')) / 7
  )) ]

  # Cut at the day of the first COVID death
  deaths_long[, in_baseline := as.integer(date_full < first_covid_death_date) ]

  # Aggregate deaths and number of days by week across all other groups
  deaths_agg <- deaths_long[,
                            .(deaths=sum(deaths), observed_days=length(unique(date_full))),
                            by=.(location_code, year, week, sex, age_group_code, in_baseline)
                            ][order(location_code, year, week, sex, age_group_code, in_baseline)]

  return(deaths_agg)
}


#' Prepare population data for Italy
#'
#' Given raw population data by age, sex, year, and province in Italy, format
#' the data for modeling. This is a convenience function specifically designed
#' to be use with IStat data. For more information, see the repository README.
#'
#' @param pop_raw Raw population data.table downloaded from IStat
#' @param age_cutoffs Vector of the starting years for each age group bin. For
#'
#' @return Data.table of formatted population data
#'
ita_prepare_pop <- function(pop_raw, age_cutoffs){
  setnames(pop_raw, c('TIME','ITTER107'), c('year', 'icode'))

  # Format sex
  pop_raw[, sex := gsub('s', '', Gender)]
  # Add province code
  pop_raw <- merge(pop_raw, location_table[, .(location_code, icode)], by='icode')
  # Format age
  pop_raw <- pop_raw[ Age != "total", ]
  pop_raw[, age_years := gsub('[ yearsandover]', '', Age) ]
  # Sort into age groups
  age_groups <- create_age_groups(age_cutoffs=age_cutoffs)
  pop_raw$age_group_code <- age_groups$age_group_code[
    findInterval(pop_raw$age_years, age_groups$age_group_years_start)
    ]

  ## Collapse by identifiers
  pop_agg <- pop_raw[,
                     .(pop = sum(Value)),
                     by=.(location_code, year, sex, age_group_code)
                     ][order(location_code, year, sex, age_group_code)]

  return(pop_agg)
}
