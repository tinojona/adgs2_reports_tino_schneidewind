#######################################################################################
#--------------------------------PHENOLOGY MODELLING----------------------------------#
#######################################################################################



# tutorial
# https://geco-bern.github.io/handfull_of_pixels/phenology_modelling.html

# pseudo mechanistic modelling:
# more hands-on than machine learning
# built in physical assumptions 
# has free parameters, that are yet undefined

# threshold for growing degree days (GDD) ?start of growing season?:
## certain number of days with average temperature above 5°C


# CHAPTER 6
# I will use the phenocamr package which 
# interfaces with the phenocam network API
# to download time series of vegetation 
# greenness and derived phenology metrics
library(phenocamr)
library(geodata)
library(terra)
library(dplyr)
library(here)
library(ggplot2)
library(patchwork)
library(GenSA)
library(BayesianTools)
library(daymetr)


# download greenness time series,
# calculate phenology (phenophases),
# amend with DAYMET data
phenocamr::download_phenocam(
  site = "harvard$",
  veg_type = "DB",
  roi_id = "1000",
  daymet = TRUE,
  phenophase = TRUE,
  trim = 2022,
  out_dir = tempdir()
)

harvard_phenocam_data <- readr::read_csv(
  file.path(tempdir(), "harvard_DB_1000_3day.csv"), 
  comment = "#"
)

# reading in harvard phenology only retaining
# spring (rising) phenology for the GCC 90th
# percentile time series (the default)
harvard_phenology <- readr::read_csv(
  file.path(
    tempdir(),
    "harvard_DB_1000_3day_transition_dates.csv"
  ),
  comment = "#"
) |>
  dplyr::filter(
    direction == "rising",
    gcc_value == "gcc_90"
  )


# return mean daily temperature as well
# as formal dates (for plotting)
harvard_temp <- harvard_phenocam_data |>
  group_by(year) |>
  dplyr::mutate(
    tmean = (tmax..deg.c. + tmin..deg.c.)/2
  ) |> 
  dplyr::mutate(
    date = as.Date(date),
    gdd = cumsum(ifelse(tmean >= 5, tmean - 5, 0))
  ) |>
  dplyr::select(
    date,
    year,
    tmean,
    gdd
  ) |>
  ungroup()

# convert the harvard phenology data and only
# retain required data
harvard_phenology <- harvard_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  select(
    year,
    doy,
    transition_25,
    threshold_25
  )





# return mean daily temperature as well
# as formal dates (for plotting)
harvard_temp <- harvard_phenocam_data |>
  group_by(year) |>
  dplyr::mutate(
    tmean = (tmax..deg.c. + tmin..deg.c.)/2
  ) |> 
  dplyr::mutate(
    date = as.Date(date),
    gdd = cumsum(ifelse(tmean >= 5, tmean - 5, 0))
  ) |>
  dplyr::select(
    date,
    year,
    tmean,
    gdd
  ) |>
  ungroup()

# convert the harvard phenology data and only
# retain required data
harvard_phenology <- harvard_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  select(
    year,
    doy,
    transition_25,
    threshold_25
  )


# growing degree of model optimisation
gdd_model <- function(temp, par) {
  # split out parameters from a simple
  # vector of parameter values
  temp_threshold <- par[1]
  gdd_crit <- par[2]
  
  # accumulate growing degree days for
  # temperature data
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0))
  
  # figure out when the number of growing
  # degree days exceeds the minimum value
  # required for leaf development, only
  # return the first value
  doy <- unlist(which(gdd >= gdd_crit)[1])
  
  return(doy)
}

# confirm that the model function
# returns expected results (i.e. DOY 114)
# (we filter out the year 2010, but
# removing the filter would run the
# model for all years!)
prediction <- harvard_temp |>
  # dplyr::filter(
  #   year == 2010
  # ) |>
  group_by(year) |>
  summarize(
    pred = gdd_model(
      temp = tmean,
      par = c(5, 130.44)
    )  
  )

print(prediction)
print(mean(prediction$pred))


# phenology model calibration
# run model and compare to true values
# returns the RMSE
rmse_gdd <- function(par, data) {
  
  # split out data
  drivers <- data$drivers
  validation <- data$validation
  
  # calculate phenology predictions
  # and put in a data frame
  predictions <- drivers |>
    group_by(year) |>
    summarise(
      predictions = gdd_model(
        temp = tmean,
        par = par
      )
    )
  
  predictions <- left_join(predictions, validation, by = "year")
  
  rmse <- predictions |>
    summarise(
      rmse = sqrt(mean((predictions - doy)^2, na.rm = TRUE))
    ) |>
    pull(rmse)
  
  # return rmse value
  return(rmse)
}


# starting model parameters
par = c(0, 130)

# limits to the parameter space
lower <- c(-10,0)
upper <- c(45,500)

# data needs to be provided in a consistent
# single data file, a nested data structure
# will therefore accept non standard data formats
data <- list(
  drivers = harvard_temp,
  validation = harvard_phenology
)

# optimize the model parameters
optim_par = GenSA::GenSA(
  par = par,
  fn = rmse_gdd,
  lower = lower,
  upper = upper,
  control = list(
    max.call = 4000
  ),
  data = data
)$par

optim_par
# After optimization routine, the optimal parameters are 2.811523, 228.27844 
# for the temperature threshold and number of accumulation days respectively



# run the model for all years
# to get the phenology predictions
predictions <- harvard_temp |>
  group_by(year) |>
  summarize(
    prediction = gdd_model(
      temp = tmean,
      par = optim_par
    )  
  )
print(predictions)



library(daymetr)

# Download daily data
# for spatial upscaling
daymetr::download_daymet_tiles(
  tiles = 11935,
  start = 2012,
  end = 2012,
  param = c('tmin','tmax'),
  path = file.path(here::here(), "/data-raw/"),
  silent = FALSE
)

# calculate the daily mean values
r <- daymetr::daymet_grid_tmean(
  path = file.path(here::here(), "/data-raw/"),
  product = 11935,
  year = 2012,
  internal = TRUE
)

# reproject to lat lon
r <- terra::project(
  r,
  "+init=epsg:4326"
)



# subset to first 180 days
ma_nh_temp <- terra::subset(
  r,
  1:180
)


predicted_phenology <- terra::app(
  ma_nh_temp,
  fun = gdd_model,
  par = optim_par
)

