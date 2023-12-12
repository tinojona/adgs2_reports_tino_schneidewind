#------------------------------------------------------
# Questions to Phenology modelling
#------------------------------------------------------
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
library(MODISTools)
library(leaflet)
library(tidyterra)


# How can you improve the model used to regionally scale the results in Chapter 6?
## provide at least three ways to improve the model

### IDEE1 optimize loss function using a different parameter! instead of rsme
### IDEE2 : add another parameter (eg radiation)
### IDEE3 : Include data from a different site to increase the models generelisability

### plus : change starting model parameters and limits to improve the computational effort of model training
### plus : include machine learning 


#########################


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


# inlcude data from another site

phenocamr::download_phenocam(
  site = "worcester",
  veg_type = "DB",
  roi_id = "1000",
  daymet = TRUE,
  phenophase = TRUE,
  trim = 2022,
  out_dir = tempdir()
)

worcester_phenocam_data <- readr::read_csv(
  file.path(tempdir(), "worcester_DB_1000_3day.csv"), 
  comment = "#"
)


worcester_phenology <- readr::read_csv(
  file.path(
    tempdir(),
    "worcester_DB_1000_3day_transition_dates.csv"
  ),
  comment = "#"
) |>
  dplyr::filter(
    direction == "rising",
    gcc_value == "gcc_90"
  )


# return mean daily temperature as well
# as formal dates (for plotting)
worcester_temp <- worcester_phenocam_data |>
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
worcester_phenology <- worcester_phenology |>
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
worcester_temp <- worcester_phenocam_data |>
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
worcester_phenology <- worcester_phenology |>
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

# Site ID einfügen
# data <- list(
#   drivers = rbind(worcester_temp, harvard_temp),
#   validation = rbind(worcester_phenology, harvard_phenology)
# )

worcester_temp$ID <- "worcester"
harvard_temp$ID <- "harvard"
worcester_phenology$ID <- "worcester"
harvard_phenology$ID <- "harvard"

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


# phenology model calibration
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
  drivers = rbind(worcester_temp, harvard_temp),
  validation = rbind(worcester_phenology, harvard_phenology)
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
# After optimization routine, the optimal parameters are 5.861427 143.048134
# for the temperature threshold and number of accumulation days respectively



# run the model for all years
# to get the phenology predictions
predictions <- rbind(worcester_temp, harvard_temp) |>
  group_by(year) |>
  summarize(
    prediction = gdd_model(
      temp = tmean,
      par = optim_par
    )  
  )
print(predictions)




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




######################
# join predicted with observed data
validation <- left_join(predictions, rbind(worcester_phenology,harvard_phenology))
# validation <- left_join(predictions, harvard_phenology)

ggplot(validation) +
  geom_smooth(
    aes(
      doy,
      prediction
    ),
    colour = "grey25",
    method = "lm"
  ) +
  geom_point(
    aes(
      doy,
      prediction
    )
  ) +
  geom_abline(
    intercept=0, 
    slope=1, 
    linetype="dotted"
  ) +
  labs(
    x = "Observed leaf-out date (DOY)",
    y = "Predicted leaf-out date (DOY)"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )


# set te colour scale manually
pal <- colorNumeric(
  "magma",
  values(predicted_phenology),
  na.color = "transparent"
)

# build the leaflet map
# using ESRI tile servers
# and the loaded demo raster
leaflet() |> 
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") |>
  addProviderTiles(providers$Esri.WorldTopoMap, group = "World Topo") |>
  addRasterImage(
    predicted_phenology,
    colors = pal,
    opacity = 0.8,
    group = "Phenology model results"
  ) |>
  addLayersControl(
    baseGroups = c("World Imagery","World Topo"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE),
    overlayGroups = c("Phenology model results")
  ) |>
  addLegend(
    pal = pal,
    values = values(predicted_phenology),
    title = "DOY")




# Statistically compare the results with the MODIS MCD12Q2 phenology product

# wir vergleichen unsere DOY "predictions" für harvard and worchester mit den Daten von MODIS
# dafür müssen wir die DOYs für unseren Standorten extrahieren, spatially mitteln,  
# und dann mittels einem statistischen test verlgleichen. (ttest wenn normalverteilt)

# download and save phenology data
phenology <- MODISTools::mt_subset(
  product = "MCD12Q2",
  lat = 42.5378,
  lon = -72.1715,
  band = "Greenup.Num_Modes_01",
  start = "2008-01-01",
  end = "2022-12-31",
  km_lr = 20,
  km_ab = 100,
  site_name = "swiss",
  internal = TRUE,
  progress = FALSE
)


# screening of data
phenology <- phenology |>
  mutate(
    value = ifelse(value > 32656, NA, value),
    value = as.numeric(format(as.Date("1970-01-01") + value, "%j")),
    value = ifelse (value < 200, value, NA)
  )
phenology_raster <- MODISTools::mt_to_terra(
  phenology,
  reproject = TRUE
)


ggplot() +
  tidyterra::geom_spatraster(data = phenology_raster) +
  scale_fill_viridis_c(
    na.value = NA,
    name = "DOY"
  ) +
  theme_bw()



phenology_filtered <- phenology |>
  select(
    value,
    calendar_date
  ) |>
  mutate(
    calendar_date = format(as.Date(calendar_date), format = "%Y")
  ) |>
  group_by(calendar_date) |>
  summarize(
    pheno_mean = mean(value, na.rm=T)
  )

plot(predictions$prediction, phenology_filtered$pheno_mean, 
     xlim = c(90,140),
     ylim = c(90, 140),
     lty = 3,
     lwd = 3,
     pch = 5,
     col = "dark green", 
     xlab = "My Predictions",
     ylab = "Modis Values")

# ZWEITER STANDORT; MITTELN VON BEIDEN; DANN VERGLEICHEN 
# statistischer vergleich
# H_0 both values are the same

t.test(predictions$prediction, phenology_filtered$pheno_mean) # p=0.0002
cor(predictions$prediction, phenology_filtered$pheno_mean) # 0.12


# join predicted with observed data
validation <- as.data.frame(cbind(predictions$prediction, phenology_filtered$pheno_mean))
colnames(validation) = c("X", "Y")

ggplot(validation) +
  geom_smooth(
    aes(
      X,
      Y
    ),
    colour = "brown",
    method = "lm"
  ) +
  geom_point(
    aes(
      X,
      Y
    )
  ) +
  geom_abline(
    intercept=0, 
    slope=1, 
    linetype="dotted"
  ) +
  labs(
    x = "My Predictions (DOY)",
    y = "MODIS (DOY)"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )


# 2013 NA
# join predicted with observed data
v <- as.data.frame(cbind(predictions$prediction, phenology_filtered$pheno_mean))
colnames(v) = c("X", "Y")
v$X[6]=NA

ggplot(v) +
  geom_smooth(
    aes(
      X,
      Y
    ),
    colour = "brown",
    method = "lm"
  ) +
  geom_point(
    aes(
      X,
      Y
    )
  ) +
  geom_abline(
    intercept=0, 
    slope=1, 
    linetype="dotted"
  ) +
  labs(
    x = "My Predictions (DOY)",
    y = "MODIS (DOY)"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )




# NOTIZEN

#
# andere site einbeziehen und daten mitteln, 
# statistische test kreativ werden




