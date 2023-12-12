#------------------------------------------------------------------------------------------------
#                LAND COVER AND CLASSIFICATION
#------------------------------------------------------------------------------------------------
rm(list=ls())
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
library(jsonlite)
library(appeears)
library(xgboost)



# set a key to the keychain
rs_set_key(
  user = "tinojona",
  password = "#Hashtag1999"
)


# token <- rs_login(user = "earth_data_user")
# 
# token <- rs_login("eyJ0eXAiOiJKV1QiLCJvcmlnaW4iOiJFYXJ0aGRhdGEgTG9naW4iLCJzaWciOiJlZGxqd3RwdWJrZXlfb3BzIiwiYWxnIjoiUlMyNTYifQ.eyJ0eXBlIjoiVXNlciIsInVpZCI6InRpbm9qb25hIiwiZXhwIjoxNzA1MTM5OTY3LCJpYXQiOjE2OTk5NTU5NjcsImlzcyI6IkVhcnRoZGF0YSBMb2dpbiJ9.HvRG2dorz_hQ8tIO0f8rgIco4E6DQT_guOxgfUu_LLVnEF-hQKe8Trq3-Myd5EOTfrx9nDwUjwTycNC46uPHMYUHZGJL9FGL7x1A9DSfyYXJgj4XRRtbEfOv0YCBudh9yW5sqT1VZkkUklyVPrtKxzRusrjaGTke5E4Cje2yEOBvn-ZtTpEWY82bhN37OXY4snVcNRciH4iAGdJD6a9sKRJ3F-oOdiqdhpJXsZu7pDl62A89p84NwLAc3HF4NrFEIj_ChY69yoRmUlGQEEB5Bb3DRta62T5lDScLaBI1eDB7QX25P43Q2oYXtG6GGeGbJRcz3X4AjL8COnh4A0zm4Q")
# 
# rs_logout(token)

# conversion from tidy data to a raster format
lai_2012 <- readRDS("./data/lai_2012.rds")

# save this data for later use
# to speed up computation
r <- MODISTools::mt_to_terra(
  lai_2012,
  reproject = TRUE
)

# convert a multi-layer raster image
# to wide dataframe
df <- as.data.frame(r, cell = TRUE)

# the content of a single feature (vector)
# limited to the first 5 values for brevity
print(df[1,1:5])

# cluster the data 
clusters <- kmeans(
  df[,-1],
  centers = 2
)



# use the original raster layout as
# a template for the new map (only
# using a single layer)
kmeans_map <- terra::rast(r, nlyr=1)

# assign to each cell value (location) of this
# new map using the previously exported cell
# values (NA values are omitted so a 1:1
# mapping would not work)
kmeans_map[df$cell] <- clusters$cluster


library(leaflet)

# set te colour scale manually
palcol <- colorFactor(
  c("#78d203", "#f9ffa4"),
  domain = 1:2,
  na.color = "transparent"
)

# build the leaflet map
leaflet() |> 
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") |>
  addProviderTiles(providers$Esri.WorldTopoMap, group = "World Topo") |>
  addRasterImage(
    kmeans_map,
    colors = palcol,
    opacity = 0.5,
    method = "ngb",
    group = "k-means cluster results"
  ) |>
  addLayersControl(
    baseGroups = c("World Imagery","World Topo"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE),
    overlayGroups = c("k-means cluster results")
  ) |>
  addLegend(
    colors = palcol(1:2),
    values = c(1, 2),
    title = "cluster",
    labels = c(1, 2)
  )

# Read the validation sites from
# Fritz et al. 2017 straight from
# Zenodo.org
validation_sites <- readr::read_csv(
  "https://zenodo.org/record/6572482/files/Global%20LULC%20reference%20data%20.csv?download=1"
)

# filter out data by competition,
# coverage percentage and latitude
# (use round brackets to enclose complex
# logical statements in a filter call!)
validation_selection <- validation_sites |>
  dplyr::filter(
    (competition == 4 | competition == 1),
    perc1 > 80,
    lat > 0
  )

# the above selection includes all data
# but we now subsample to 150 random locations
# per (group_by()) land cover class (LC1)
# set a seed for reproducibilty
set.seed(0)

validation_selection <- validation_selection |>
  dplyr::slice_sample(n = 150, by = LC1)

# split validation selection
# by land cover type into a nested
# list, for easier processing
# later on
validation_selection <- validation_selection |>
  dplyr::group_by(LC1) |>
  dplyr::group_split()


#!!!! 
# for every row download the data for this
# location and the specified reflectance
# bands
task_nbar <- lapply(validation_selection, function(x){
  
  # loop over all list items (i.e. land cover classes)
  base_query <- x |>
    dplyr::rowwise() |>
    do({
      data.frame(
        task = paste0("nbar_lc_",.$LC1),
        subtask = as.character(.$pixelID),
        latitude = .$lat,
        longitude = .$lon,
        start = "2012-01-01",
        end = "2012-12-31",
        product = "MCD43A4.061",
        layer = paste0("Nadir_Reflectance_Band", 1:4)
      )
    }) |>
    dplyr::ungroup()
  
  # build a task JSON string 
  task <- rs_build_task(
    df = base_query
  )
  
  # return task
  return(task)
})

# Query the appeears API and process
# data in batches - this function
# requires an active API session/login
rs_request_batch(
  request = task_nbar,
  workers = 10,
  user = "tinojona",
  path = tempdir(),
  verbose = TRUE,
  time_out = 28800
)


# list all MCD43A4 files, note that
# that list.files() uses regular
# expressions when using wildcards
# such as *, you can convert general
# wildcard use to regular expressions
# with glob2rx()
files <- list.files(
  tempdir(),
  glob2rx("*MCD43A4-061-results*"),
  recursive = TRUE,
  full.names = TRUE
)

# read in the data (very fast)
# with {vroom} and set all
# fill values (>=32767) to NA
nbar <- vroom::vroom(files)
nbar[nbar >= 32767] <- NA

# retain the required data only
# and convert to a wide format
nbar_wide <- nbar |>
  dplyr::select(
    Category,
    ID,
    Date,
    Latitude,
    Longitude,
    starts_with("MCD43A4_061_Nadir")
  ) |>
  tidyr::pivot_wider(
    values_from = starts_with("MCD43A4_061_Nadir"),
    names_from = Date
  )

# split out only the site name,
# and land cover class from the
# selection of validation sites
# (this is a nested list so we
# bind_rows across the list)
sites <- validation_selection |>
  dplyr::bind_rows() |>
  dplyr::select(
    pixelID,
    LC1
  ) |>
  dplyr::rename(
    Category = "pixelID"
  )

# combine the NBAR and land-use
# land-cover labels by location
# id (Category)
ml_df <- left_join(nbar_wide, sites) |>
  dplyr::select(
    LC1,
    contains("band")
  )

# select packages
# avoiding tidy catch alls
library(rsample)

# create a data split across
# land cover classes
ml_df_split <- ml_df |>
  rsample::initial_split(
    strata = LC1,
    prop = 0.8
  )

# select training and testing
# data based on this split
train <- readRDS("./data/training_data.rds")
test <- readRDS("./data/test_data.rds")

# train <- rsample::training(ml_df_split)
# test <- rsample::testing(ml_df_split)


# load the parsnip package
# for tidy machine learning
# modelling and workflows
# to manage workflows
library(parsnip)
library(workflows)

# specify our model structure
# the model to be used and the
# type of task we want to evaluate
model_settings <- parsnip::boost_tree(
  trees = 50,
  min_n = tune(),
  tree_depth = tune(),
  # learn_rate = tune()
) |>
  set_engine("xgboost") |>
  set_mode("classification")

# create a workflow compatible with
# the {tune} package which combines
# model settings with the desired
# model structure (data / formula)
xgb_workflow <- workflows::workflow() |>
  add_formula(as.factor(LC1) ~ .) |>
  add_model(model_settings)

print(xgb_workflow)





# load the dials package
# responsible for (hyper) parameter
# sampling schemes to tune
# parameters (as extracted)
# from the model specifications
library(tune)
library(dials)

hp_settings <- dials::grid_latin_hypercube(
  tune::extract_parameter_set_dials(xgb_workflow),
  size = 3
)

print(hp_settings)



# set the folds (division into different)
# cross-validation training datasets
folds <- rsample::vfold_cv(train, v = 3)

# optimize the model (hyper) parameters
# using the:
# 1. workflow (i.e. model)
# 2. the cross-validation across training data
# 3. the (hyper) parameter specifications
# all data are saved for evaluation
xgb_results <- tune::tune_grid(
  xgb_workflow,
  resamples = folds,
  grid = hp_settings,
  control = tune::control_grid(save_pred = TRUE)
)


# select the best model based upon
# the root mean squared error
xgb_best <- tune::select_best(
  xgb_results,
  metric = "roc_auc"
)

# cook up a model using finalize_workflow
# which takes workflow (model) specifications
# and combines it with optimal model
# parameters into a model workflow
xgb_best_hp <- tune::finalize_workflow(
  xgb_workflow,
  xgb_best
)

print(xgb_best_hp)

# train a final (best) model with optimal
# hyper-parameters
xgb_best_model <- fit(xgb_best_hp, train)



# run the model on our test data
# using predict()
test_results <- predict(xgb_best_model, test)

# load the caret library to
# access confusionMatrix functionality
library(caret)
library(lattice)
# use caret's confusionMatrix function to get
# a full overview of metrics
caret::confusionMatrix(
  reference = as.factor(test$LC1),
  data = as.factor(test_results$.pred_class)
)



# We can define an appeears
# download task using a simple
# dataframe and a map from which
# an extent is extracted
task_df <- data.frame(
  task = "raster_download",
  subtask = "swiss",
  start = "2012-01-01",
  end = "2012-12-31",
  product = "MCD43A4.061",
  layer = paste0("Nadir_Reflectance_Band", 1:4)
)

# build the area based request/task
# using the extent of our previous
# kmeans map, export all results
# as geotiff (rather than netcdf)
task <- rs_build_task(
  df = task_df,
  roi = kmeans_map,
  format = "geotiff"
)

# request the task to be executed
# with results stored in a
# temporary location (can be changed)
rs_request(
  request = task,
  user = "your_api_id",
  transfer = TRUE,
  path = tempdir(),
  verbose = TRUE
)


files <- list.files(
  tempdir(),
  "*Reflectance*",
  recursive = TRUE,
  full.names = TRUE
)

# load this spatial data to run the model
# spatially
swiss_multispec_data <- terra::rast(files)


# the model only works when variable names
# are consistent we therefore rename them
band_names <- data.frame(
  name = names(swiss_multispec_data)
) |>
  mutate(
    date = as.Date(substr(name, 40, 46), format = "%Y%j"),
    name = paste(substr(name, 1, 35), date, sep = "_"),
    name = gsub("\\.","_", name)
  )

# reassign the names of the terra image stack
names(swiss_multispec_data) <- band_names$name



# return probabilities, where each class is
# associated with a layer in an image stack
# and the probabilities reflect the probabilities
# of the classification for said layer
lulc_probabilities <- terra::predict(
  swiss_multispec_data,
  xgb_best_model,
  type = "prob"
)


# generate the map by selecting maximum probabilities
# from the model output
lulc_map <- terra::app(lulc_probabilities, which.max)






classes <- c(
  "Tree Cover",
  "Shrub Cover",
  "Herbaceous Vegetation & Grassland",
  "Cultivated and Managed",
  "Mosaic: Managed & Natural Vegetation",
  "Regularly Flooded & Wetland",
  "Urban & Built Up",
  "Snow and Ice",
  "Barren",
  "Open Water"
)

# set te colour scale manually
palcol <- colorFactor(
  c(
    "#05450a",
    "#78d203",
    "#009900",
    "#c24f44",
    "#ff6d4c",
    "#27ff87",
    "#a5a5a5",
    "#69fff8",
    "#f9ffa4",
    "#1c0dff"
  ),
  na.color = NA,
  domain = 1:10
)

# build the leaflet map
leaflet() |> 
  addProviderTiles(providers$Esri.WorldImagery, group = "World Imagery") |>
  addProviderTiles(providers$Esri.WorldTopoMap, group = "World Topo") |>
  addRasterImage(
    lulc_map,
    colors = palcol,
    opacity = 0.8,
    method = "ngb",
    group = "XGBOOST"
  ) |>
  addRasterImage(
    modis_lulc,
    colors = palcol,
    opacity = 0.8,
    method = "ngb",
    group = "MODIS MCD12Q1"
  ) |>
  addLayersControl(
    baseGroups = c("World Imagery","World Topo"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE),
    overlayGroups = c("XGBOOST", "MODIS MCD12Q1")
  ) |>
  hideGroup("MODIS MCD12Q1") |>
  addLegend(
    colors = palcol(1:10),
    values = 1:10,
    labels = classes,
    title = "Land-Use and Land-Cover class"
  )