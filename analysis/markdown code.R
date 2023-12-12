

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
library(raster)
library(MODISTools)
library(leaflet)
library(tidyterra)



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

harvard_phenology <- harvard_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  dplyr::select(
    year,
    doy,
    transition_25,
    threshold_25
  )

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

harvard_phenology <- harvard_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  dplyr::select(
    year,
    doy,
    transition_25,
    threshold_25
  )







source("./functions/phenology_models.R")








# starting model parameters
par = c(0, 130)
lower <- c(-10,0)
upper <- c(45,500)

data_harvard <- list(
  drivers = harvard_temp,
  validation = harvard_phenology
)

optim_par_harvard = GenSA::GenSA(
  par = par,
  fn = rmse_gdd,
  lower = lower,
  upper = upper,
  control = list(
    max.call = 4000
  ),
  data = data_harvard
)$par








predictions_harvard <- harvard_temp |>
  group_by(year) |>
  summarize(
    prediction = gdd_model(
      temp = tmean,
      par = optim_par_harvard
    )  
  )




#hier


prediction_harvard_10 <- harvard_temp |>
  dplyr::filter(
    year == 2010
  ) |>
  group_by(year) |>
  summarize(
    pred = gdd_model(
      temp = tmean,
      par = optim_par_harvard
    )  
  )

# Download daily data
daymetr::download_daymet_tiles(
  tiles = 11935,
  start = 2010,
  end = 2010,
  param = c("tmin","tmax"),
  path = paste0(here::here(), "/data-raw/"),
  silent = TRUE
)

# calculate the daily mean values
r <- daymetr::daymet_grid_tmean(
  path = paste0(here::here(), "/data-raw/"),
  product = 11935,
  year = 2010,
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

predicted_phenology_harvard <- terra::app(
  ma_nh_temp,
  fun = gdd_model,
  par = optim_par_harvard
)










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

worcester_phenology <- worcester_phenology %>%
  dplyr::mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) # hier war eigentlich |>
worcester_phenology <- readr::read_csv(
  file.path(
    tempdir(),
    "worcester_DB_1000_3day_transition_dates.csv"
  ),
  comment = "#"
) %>%
  dplyr::filter(
    direction == "rising",
    gcc_value == "gcc_90"
  )

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

worcester_phenology <- worcester_phenology |>
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")),
    year = as.numeric(format(as.Date(transition_25),"%Y"))
  ) |>
  dplyr::select(
    year,
    doy,
    transition_25,
    threshold_25
  )



# Site ID
worcester_temp$ID <- "worcester"
harvard_temp$ID <- "harvard"
worcester_phenology$ID <- "worcester"
harvard_phenology$ID <- "harvard"
#hier


data_ha_wo <- list(
  drivers = rbind(worcester_temp, harvard_temp),
  validation = rbind(worcester_phenology, harvard_phenology)
)

optim_par_ha_wo = GenSA::GenSA(
  par = par,
  fn = rmse_gdd,
  lower = lower,
  upper = upper,
  control = list(
    max.call = 4000
  ),
  data = data_ha_wo
)$par

predictions_ha_wo <- rbind(worcester_temp, harvard_temp) |>
  group_by(year) |>
  summarize(
    prediction = gdd_model(
      temp = tmean,
      par = optim_par_ha_wo
    )  
  )





prediction_ha_wo_10 <- rbind(worcester_temp, harvard_temp) |>
  dplyr::filter(
    year == 2010
  ) |>
  group_by(year) |>
  summarize(
    pred = gdd_model(
      temp = tmean,
      par = optim_par_ha_wo
    )  
  )


predicted_phenology_ha_wo <- terra::app(
  ma_nh_temp,
  fun = gdd_model,
  par = optim_par_ha_wo
)





phenology <- MODISTools::mt_subset(
  product = "MCD12Q2",
  lat = 42.5378,
  lon = -72.1715,
  band = "Greenup.Num_Modes_01",
  start = "2008-01-01",
  end = "2022-12-31",
  km_lr = 20,
  km_ab = 20,
  site_name = "swiss",
  internal = TRUE,
  progress = FALSE
)






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

phenology_filtered <- phenology |>
  dplyr::select(
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





phenology <- MODISTools::mt_subset(
  product = "MCD12Q2",
  lat = 42.2697,
  lon = -71.8428,
  band = "Greenup.Num_Modes_01",
  start = "2008-01-01",
  end = "2022-12-31",
  km_lr = 20,
  km_ab = 20,
  site_name = "swiss",
  internal = TRUE,
  progress = FALSE
)











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

phenology_filtered_2 <- phenology |>
  dplyr::select(
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

phenology_filtered$mean_2 <- rep(1, length(phenology_filtered$calendar_date))

for (i in 1:length(phenology_filtered$calendar_date)){
  phenology_filtered$mean_2[i] = mean(c(phenology_filtered$pheno_mean[i], phenology_filtered_2$pheno_mean[i]))
}








phenology <- MODISTools::mt_subset(
  product = "MCD12Q2",
  lat = 42.5378,
  lon = -72.1715,
  band = "Greenup.Num_Modes_01",
  start = "2010-01-01",
  end = "2010-12-31",
  km_lr = 20,
  km_ab = 20,
  site_name = "swiss",
  internal = TRUE,
  progress = FALSE
)
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




##########################################################################
# compare two rasters




# predicted_phenology_ha_wo (my raster)
# phenology_raster (MODIS raster)
# raster_pred_modis is my predicted phenology ha wo croped to the size of modis raster

raster_pred_modis <- crop(predicted_phenology_ha_wo,phenology_raster)

# ggplot() +
#   tidyterra::geom_spatraster(data = predicted_phenology_ha_wo) +
#   scale_fill_viridis_c(
#     na.value = NA,
#     name = "DOY"
#     ) +
#   theme_bw()
# 
# 
# ggplot() +
#   tidyterra::geom_spatraster(data = raster_pred_modis) +
#   scale_fill_viridis_c(
#     na.value = NA,
#     name = "DOY"
#     ) +
#   theme_bw()


pred_resampled <- terra::resample(raster_pred_modis,phenology_raster)



plot(pred_resampled,phenology_raster,
     lwd=4, 
     xlab = "MODIS (DOY)",
     ylab = "My Predictions (DOY)")
     #xlim = c(80,110),
     #ylim = c(110, 125))
abline(a=79.40996, b=0.39753, col="brown1")


raster_predictions <- raster::as.matrix(pred_resampled)
raster_modis       <- raster::as.matrix(phenology_raster)
linmod <- lm(raster_predictions ~ raster_modis)

summary(linmod)




validation <- as.data.frame(cbind(raster_modis, raster_predictions))

colnames(validation) = c("X", "Y")

ggplot(validation) +
  geom_smooth(
    aes(
      Y,
      X
    ),
    colour = "brown2",
    method = "lm"
  ) +
  geom_point(
    aes(
      Y,
      X
    )
  ) +
  geom_abline(
    intercept=0, 
    slope=1, 
    linetype="dotted"
  ) +
  labs(
    x = "My Predictions (DOY)",
    y = "MODIS (DOY)",
    title = "Agreement between predictions and MODIS observations for harvard"
  ) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )




# 
library(diffeR)

categoryComponentsPlot( comp= pred_resampled , ref = phenology_raster )
