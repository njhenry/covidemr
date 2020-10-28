
#' Create age group categories
#'
#' @description Given a set of age cutoffs starting at 0, this function
#'   deterministically creates a data.table of age groupings with IDs and
#'   lower/upper bounds for each age group
#'
#' @param age_cutoffs Vector of lower bounds for each age group, starting at 0
#'
#' @return Data.table containing age group categories
#'
#' @import data.table
#' @export
create_age_groups <- function(age_cutoffs = seq(0, 100, by=20)){
  if(length(age_cutoffs)==0) stop("Need at least one age group!")
  age_groups <- data.table::data.table(age_group_code = 1:length(age_cutoffs))
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


#' Read CSV downloaded from iStat
#'
#' @description Read a CSV downloaded from iStat, the Italian Statistical
#'   Agency, as an R data.table
#'
#' @param fp [char] The filepath of the input file to read
#'
#' @return A data.table representation of the CSV, with NAs appropriately placed
#'
#' @import data.table
#' @export
ita_read_istat_csv <- function(fp){
  data.table::fread(fp, na.strings=c("","NA","n.d."))
}

#' Prepare a single year of death data for Italy
#'
#' @description Prepare a single year of death data for Italy. This function is
#'   wrapped by `ita_prepare_deaths_all_years()`.
#'
#' @param deaths_raw Data.table downloaded from IStat
#' @param age_cutoffs Vector of the starting years for each age group bin. For
#'   more information, see \link{\code{create_age_groups}}
#' @param this_year Year to be extracted
#' @param first_covid_death_date R Date object designating the date of the first
#'   confirmed COVID death in Italy.
#'
#' @return Data.table of formatted death data
#'
#' @import data.table
#' @export
ita_prepare_deaths_single_year <- function(
  deaths_raw, age_cutoffs, this_year, first_covid_death_date
){
  # Check data validity
  if(length(this_year) != 1) stop(
    "Parameter this_year should take a single value. To prepare multiple years",
    "of data, use the wrapper function ita_prepare_deaths_all_years()."
  )

  # Keep only required columns and reshape long
  id_cols <- c('PROV','CL_ETA','GE')
  val_cols <- paste0(c("M_", "F_"), this_year - 2000)
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
  deaths_long[, year := this_year ]
  # Add month and day, then convert to week ID
  deaths_long[, month := as.integer(substr(day_id, 1, 2))]
  deaths_long[, day := as.integer(substr(day_id, 3, 4)) ]
  deaths_long[, date_full := as.Date(
    sprintf('%04d-%02d-%02d', year, month, day),
    format = '%Y-%m-%d'
  )]
  # Weeks are defined as starting on January 1st (Jan 8 starts week 2, etc.)
  # Set days 365 and 366 to week 52 to avoid a stub week
  deaths_long[, week := as.integer(ceiling(as.numeric(strftime(date_full, format='%j') ) / 7))]
  deaths_long[ week > 52, week := 52]

  # Cut at the day of the first COVID death
  deaths_long[, in_baseline := as.integer(date_full < first_covid_death_date) ]

  # Aggregate deaths across all other groups
  deaths_agg <- deaths_long[
    , .(deaths=sum(deaths)),
    by=.(location_code, year, week, sex, age_group_code, in_baseline)
  ][order(location_code, year, week, sex, age_group_code, in_baseline)]

  # Get days in week for each grouping
  days_in_week <- data.table(date_full = seq(
    as.Date(paste0(this_year,'-01-01')),
    as.Date(paste0(this_year,'-12-31')),
    by='1 day'
  ))
  days_in_week <- days_in_week[ date_full <= max(na.omit(deaths_long$date_full)), ]
  days_in_week[, year := this_year ]
  days_in_week[, week := as.integer(ceiling(as.numeric(strftime(date_full, format='%j') ) / 7))]
  days_in_week[ week > 52, week := 52]
  days_in_week[, in_baseline := as.integer(date_full < first_covid_death_date) ]
  days_in_week_agg <- days_in_week[, .(observed_days = .N), by=.(year, week, in_baseline)]

  # Merge onto aggregated dataset
  deaths_agg <- merge(
    deaths_agg,
    days_in_week_agg,
    by = c('year','week','in_baseline')
  )

  return(deaths_agg)
}


#' Prepare death data for Italy
#'
#' @description Given raw data all-cause mortality by age, sex, date, and
#'   province in Italy, prepare the data using standard identifiers and grouped
#'   by year and week. This function is designed specifically to work with data
#'   downloaded from Istat. This is a wrapper for the function
#'   `ita_prepare_deaths_single_year()`. For more information, see the
#'   repository README.
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
#' @import data.table
#' @export
ita_prepare_deaths_all_years <- function(
  deaths_raw, age_cutoffs, model_years, first_covid_death_date
){
  # Iteratively prepare deaths to reduce memory use
  deaths_prepped_list <- vector('list', length = length(model_years))

  for(year_idx in 1:length(model_years)){
    this_year <- model_years[year_idx]
    message("Preparing deaths for year ", this_year)
    deaths_prepped_list[[year_idx]] <- ita_prepare_deaths_single_year(
      deaths_raw = deaths_raw,
      age_cutoffs = age_cutoffs,
      this_year = this_year,
      first_covid_death_date = first_covid_death_date
    )
  }

  deaths_prepped <- rbindlist(deaths_prepped_list)
  return(deaths_prepped)
 }


#' Prepare population data for Italy
#'
#' @description Given raw population data by age, sex, year, and province in
#'   Italy, format the data for modeling. This is a convenience function
#'   specifically designed to be use with IStat data. For more information, see
#'   the repository README.
#'
#' @param pop_raw Raw population data.table downloaded from IStat
#' @param age_cutoffs Vector of the starting years for each age group bin. For
#'
#' @return Data.table of formatted population data
#'
#' @import data.table
#' @export
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
  # Apply fix for Sud Sardegna province
  pop_raw <- backfill_input_data(pop_raw, index_field='icode', check_vals='IT111')

  # Collapse by identifiers
  pop_agg <- pop_raw[,
                     .(pop = sum(Value)),
                     by=.(location_code, year, sex, age_group_code)
                     ][order(location_code, year, sex, age_group_code)]

  return(pop_agg)
}


#' Backfill missing historical data
#'
#' @description Fill missing historical estimates of data using the earliest
#'   available data value.
#'
#' @details This fix is relevant for administrative boundaries that have changed
#'   over time. At the province level, this is relevant for Sud Sardegna
#'   province, which is generally listed as its previous component provinces
#'   through 2016. Sud Sardegna is a composite of territories from several
#'   now-defunct provinces, including the full territories of Medio Campidano
#'   and Carbonia-Iglesias provinces. Because there is not an exact
#'   correspondence with provinces before 2017, replicate the oldest available
#'   data for Sud Sardegna to previous years.
#'
#'   The same issue applies to data preparation for the communes (comuna), which
#'   have less stable boundaries than regions and provinces. Currently in the
#'   data preparation workflow, preparing population is required to match the
#'   subset of communes available for deaths.
#'
#' @param input_data Input data.table containing input data. The field `year`
#'   (int) is expected
#' @param index_field Field containing unique geographic identifiers by location
#' @param check_vals Values of the index field to check for backfilling
#'
#' @return Updated data.table with values for Sud Sardegna across all years
#'
#' @import data.table
#' @export
backfill_input_data <- function(input_data, index_field, check_vals){
  # Determine which years the fix might be needed for
  all_years <- sort(unique(input_data$year))

  # Create list to append backfilled values to
  concat_list <- lapply(1:(length(all_years) * length(check_vals)), function(x) NULL)
  concat_list[[1]] <- copy(input_data)
  fill_idx <- 2

  # Fill values as needed
  for(check_val in check_vals){
    if(!check_val %in% input_data[[index_field]]){
      stop(paste("Missing index value", check_val))
    }
    first_data_year <- input_data[ get(index_field) == check_val, min(year) ]
    backfill_years <- all_years[all_years < first_data_year]
    if(length(backfill_years) > 0){
      message(
        "    Filling ", check_val, "for years: ", paste(backfill_years, collapse=', ')
      )
    }
    for(backfill_year in backfill_years){
      # Copy earliest year of data and change the year to the backfill year
      concat_list[[fill_idx]] <- copy(
        input_data[(get(index_field) == check_val) & (year==first_data_year),]
      )
      concat_list[[fill_idx]][, year := backfill_year ]
      fill_idx <- fill_idx + 1
    }
  }

  # Combine into the final data.table and return
  filled_data <- rbindlist(concat_list)[order(year)]
  return(filled_data)
}


#' Create template data.table to store all observed data
#'
#' @description Create a data table containing all location-age-sex-year-week
#'   groups that will be estimated in the model. This template can then be
#'   left-joined to death, population, and covariate data.
#'
#' @param model_years Vector of years to be modeled
#' @param model_week_range Vector of weeks IDs (1-52) to include in the model
#' @param location_table Table with `location_code`s for all locations to be
#'   modeled
#' @param age_groups Age group index, created by `create_age_groups()`
#'
#' @return Template data.table containing all combinations of these indexes
#'
#' @import data.table
#' @export
create_template_dt <- function(
  model_years, model_week_range, location_table, age_groups
){
  # Create indices for age, sex, year, and week
  sex_index <- data.table(sex = c('male','female'))[, idx_sex := .I - 1 ]
  years_index <- data.table(year = sort(model_years))[, idx_year := .I - 1 ]
  week_dt <- data.table(
    week = min(model_week_range):max(model_week_range)
  )[, idx_week := .I - 1 ]
  age_groups[, idx_age := .I - 1 ]

  # Full template dataset contains sex, location, year, week, age
  # idx_<index> is a zero-indexed code for each ID column
  template_dt <- location_table[order(location_code)
    ][, `:=` (c_j = 1, idx_loc = .I - 1)
    ][ sex_index[, c_j := 1], on='c_j', allow.cartesian=TRUE
    ][ age_groups[, c_j := 1], on='c_j', allow.cartesian=TRUE
    ][ years_index[, c_j := 1], on='c_j', allow.cartesian=TRUE
    ][ week_dt[, c_j := 1], on='c_j', allow.cartesian=TRUE
    ][, c_j := NULL ]

  return(template_dt)
}
