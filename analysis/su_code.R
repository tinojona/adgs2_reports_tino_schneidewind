############################################################################################
## SPATIAL UPSCALING CODE 
## we test the reliability of 
## models predicting leaf nitrogen

# for help see ADGS1 Chapter 10+11

# 0 LIBRARIES ----------------------------
rm(list=ls())


library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(raster)
library(leaflet)
library(tidyterra)
library(ranger)
library(lattice)
library(caret)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(hrbrthemes)

# 0 Functions -------------------

# make model evaluation into a function to reuse code
eval_model <- function(mod, df_train, df_test){
  
  # add predictions to the data frames
  df_train <- df_train |> 
    drop_na()
  df_train$fitted <- predict(mod, newdata = df_train)
  
  df_test <- df_test |> 
    drop_na()
  df_test$fitted <- predict(mod, newdata = df_test)
  
  # get metrics tables
  metrics_train <- df_train |> 
    yardstick::metrics(leafN, fitted)
  
  metrics_test <- df_test |> 
    yardstick::metrics(leafN, fitted)
  
  # extract values from metrics tables
  rmse_train <- metrics_train |> 
    filter(.metric == "rmse") |> 
    pull(.estimate)
  rsq_train <- metrics_train |> 
    filter(.metric == "rsq") |> 
    pull(.estimate)
  
  rmse_test <- metrics_test |> 
    filter(.metric == "rmse") |> 
    pull(.estimate)
  rsq_test <- metrics_test |> 
    filter(.metric == "rsq") |> 
    pull(.estimate)
  
  # visualise as a scatterplot
  # adding information of metrics as sub-titles
  plot_1 <- ggplot(data = df_train, aes(leafN, fitted)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    labs(subtitle = bquote( italic(R)^2 == .(format(rsq_train, digits = 2)) ~~
                              RMSE == .(format(rmse_train, digits = 3))),
         title = "Training set") +
    theme_classic()
  
  plot_2 <- ggplot(data = df_test, aes(leafN, fitted)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    labs(subtitle = bquote( italic(R)^2 == .(format(rsq_test, digits = 2)) ~~
                              RMSE == .(format(rmse_test, digits = 3))),
         title = "Test set") +
    theme_classic()
  
  out <- cowplot::plot_grid(plot_1, plot_2)
  
  return(out)
}

# 1 THE DATA -------------------------------
## leaf N content is vital to understanding 
## photosynthesis rates and biogeochem cycles
df <- readr::read_csv("https://raw.githubusercontent.com/stineb/leafnp_data/main/data/leafnp_tian_et_al.csv")

## We will work with a limited subset of the 
## variables available in the file, and with 
## the data aggregated by sites (identified
## by their respective longitudes and latitudes):
common_species <- df |> # get species names
  group_by(Species) |> 
  summarise(count = n()) |> 
  arrange(desc(count)) |> 
  slice(1:50) |> 
  pull(Species)

dfs <- df |>   # select only needed variables
  dplyr::select(leafN, lon, lat, elv, mat, map, ndep, mai, Species) |> 
  filter(Species %in% common_species)
# group_by(lon, lat) |> 
# summarise(across(where(is.numeric), mean))

# quick overview of data
skimr::skim(dfs)

# 2.1 EXERCISES ------------------------------
## Read the paper by Ludwig et al. 2023

## 1) 
## Explain the difference between a random cross-validation 
## and a spatial cross-validation.

### There are differences in the variable selection in both methods.
### Random cross validation is able to reproduce clustered data well,
### but might not be suitable for spatial scaling due to:
### 1) a strong clustering of the reference data 
### 2) issues with quality control for the spatial predictions
### 3) Therefore it hard to quality check predictions at a new sample site
### The advantage of spatial cross-validation is, that
### 1) it works better with clustered reference data because
### it focuses on choosing variables that improve the 
### predictability of new areas.


## 2) 
## In spatial upscaling, we model
## based on environmental covariates. This implies
## that we assume the training data to sufficiently
## represent the conditions on which the model
## will be applied for generating predictions.
## Prediction errors may increase with an increasing
## distance of the prediction location from
## the training locations. The paper by Ludwig(2023)
## considers this “distance” as a geographical distance
## in Euclidian space. 
## Do you see an alternative to measuring a distance
## that considers the task of spatial upscaling
## based on environmental covariates more directly?

### We could use the distance between environmental covariates respective to their means
### (statistical distance)


# 2.2 RANDOM CROSS VALIDATION ---------------------------------------

# model formulation
mf <- recipes::recipe(leafN ~ elv + mat + map + ndep + mai + Species, 
                      data = dfs) |> 
  recipes::step_center(recipes::all_numeric(), -recipes::all_outcomes()) |>
  recipes::step_scale(recipes::all_numeric(), -recipes::all_outcomes())

# random forest model
mod <- train(
  mf, 
  data = dfs %>%
    drop_na(), 
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
  tuneGrid = expand.grid( .mtry = 3,
                          .min.node.size = 12,
                          .splitrule = "variance"),
  metric = "RMSE",
  replace = FALSE,
  sample.fraction = 0.5,
  num.trees = 2000,           # high number ok since no hperparam tuning
  seed = 1982                # for reproducibility
)


# see model performance, R2 and RSME
print(mod$results[,4:5])


# get coast outline
coast <- rnaturalearth::ne_coastline(scale = 110, returnclass = "sf")

# map of leafN distribution
ggplot() +
  
  # plot coastline
  geom_sf(data = coast,
          colour = 'black',
          size = 0.2) +
  
  # set extent in longitude and latitude
  coord_sf(
    ylim = c(-60, 80),
    expand = FALSE) +  # to draw map strictly bounded by the specified extent
  
  # plot points on map
  geom_point(data = dfs, aes(x = lon, y = lat), color = "red", size = 0.4) +
  labs(x = "", y = "") +
  theme(legend.position = "bottom")




# 2.3 SPATIAL CROSS VALIDATION ----------------------------------------

## 1) 
## What do you observe? 
## Discuss the potential implications of the geographical distribution 
## of data points for spatial upscaling.

### We can observe severe clustering of sample sites in Europe and China. 
### We only have very few samples from the US (for NA) and northern SA.
### And only 2 sample sites in Australia and 1 in Afrika. 

### This implicates that spatial cross validation is advantageous compared to 
### random cross validation. It will perform better when describing Afrika, Australia and Russia.,
### because it will chose variables which focus on a higher quality output for new sites. 


## 2) 
## Perform a spatial cross-validation. To do so, 
## first identify geographical clusters of the data using the k-means algorithm 
## (an unsupervised machine learning method), 
## considering the longitude and latitude of data points and setting k = 5.
## Plot points on a global map, showing the five clusters with distinct colors.

# cluster the data 
set.seed(1234)
dfs$clust <- as.factor(kmeans(
  dfs[,2:3],
  centers = 5
)$cluster)

# colors for plot
my_colors = c("brown3", "gold", "green2", "skyblue3","purple")

# Map of leafN ditribution with clusters
ggplot() +
  # plot coastline
  geom_sf(data = coast,
          colour = 'black',
          size = 0.2) +
  
  # set extent in longitude and latitude
  coord_sf(
    ylim = c(-60, 80),
    expand = FALSE) +  # to draw map strictly bounded by the specified extent
  
  # plot points on map
  geom_point(data = dfs, aes(x = lon, y = lat, color= clust), size = 0.6) +
  scale_color_manual(values=my_colors) +
  labs(x = "", y = "") +
  theme(legend.position = "bottom") +
  labs(title="Modelinput: lon+lat")




## 3) 
## Plot the distribution of leaf N by cluster.
dfs |>
  ggplot( aes(x=clust, y=leafN, fill=clust)) +
  geom_boxplot() +
  scale_fill_manual(values=my_colors) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")



dfs |>
  ggplot(aes(x=leafN)) +
  geom_histogram(fill = "brown3", color = "black") +
  facet_wrap(~clust) 
  

## 4) 
## Split your data into five folds that correspond to the geographical clusters 
## identified by in (2.), and fit a random forest model with the same hyperparameters 
## as above and performing a 5-fold cross-validation with the clusters as folds. 
## Report the RMSE and the R2 determined on each of the five folds

# create folds based on clusters
# assuming 'df' contains the data and a column called 'cluster' containing the 
# result of the k-means clustering
group_folds_train <- purrr::map(
  seq(length(unique(dfs$clust))),
  ~ {
    dfs |> 
      select(clust) |> 
      mutate(idx = 1:n()) |> 
      filter(clust != .) |> 
      pull(idx)
  }
)

group_folds_test <- purrr::map(
  seq(length(unique(dfs$clust))),
  ~ {
    dfs |> 
      select(clust) |> 
      mutate(idx = 1:n()) |> 
      filter(clust == .) |> 
      pull(idx)
  }
)




# create a function that trains a random forest model on a given set of rows and 
# predicts on a disjunct set of rows
train_test_by_fold <- function(df, idx_train, idx_val){ # df gelöscht
  
  mod <- ranger::ranger(
    x =  df[idx_train,2:9],  # data frame with columns corresponding to predictors, mit eckige klammern ist richtig
    y =  df$leafN[idx_train]   # a vector of the target values (not a data frame!)
  )
  
  pred <- predict(mod,       # the fitted model object 
                  data = df[idx_val,2:9] # a data frame with columns corresponding to predictors
  )
  
  # df$pred <- pred$predictions # diese zeile war eigentlich nicht da
  
  rsq <-  summary(lm(pred$predictions~df$leafN[idx_val]))$r.squared   
              # the R-squared determined on the validation set
  
  rmse <- sqrt( mean (( df$leafN[idx_val]-pred$predictions )^2 )) 
              # the root mean square error on the validation set
  
  return(tibble(rsq = rsq, rmse = rmse))
}

# apply function on each custom fold and collect validation results in a nice
# data frame
out <- purrr::map2_dfr(
  group_folds_train,
  group_folds_test,
  ~train_test_by_fold(dfs,.x, .y) # das muss man hier definieren
) |> 
  mutate(test_fold = 1:5)
out


# 
# 
# for (i in 1:5) {
#   out <- purrr::map2_dfr(
#     group_folds_train,
#     group_folds_test,
#     ~train_test_by_fold(dfs, group_folds_train[[i]],group_folds_test[[i]]) # das muss man hier definieren
#   ) #|> 
#     #mutate(test_fold = 1:5)
#   print(out)
# }
#  
# 
# unlist(group_folds_test[[3]])
# 
# mod <- ranger::ranger(
#   x =  dfs[,2:9],  # data frame with columns corresponding to predictors, mit eckige klammern ist richtig
#   y =  dfs$leafN)
# 
#   
#   pred <- predict(mod,       # the fitted model object 
#                   data = dfs[,2:9] # a data frame with columns corresponding to predictors
#   )
#   
# pred$predictions
#   rsq <-  summary(lm(pred$predictions~dfs$leafN))$r.squared   
#   # the R-squared determined on the validation set
#   
#   rmse <- sqrt(mean(lm(unname(pred)~idx_val$leafN)$residuals^2)) 
#   # the root mean square error on the validation set
#   
#   
  
  
  
  
  
  
  
  
## 5) Compare the results of the spatial cross-validation to the results of the 
## random cross-validation and discuss reasons for why you observe a difference in the 
## cross-validation metrics (if you do).







# 2.4 ENVIRONMENTAL CROSS VALIDATION -----------------------------------

## 1) To do so, perform a custom cross-validation as above, but this time 
## considering five clusters of points not in geographical space, 
## but in environmental space - spanned by the mean annual precipitation and 
## the mean annual temperature. Report the R-squared and the RMSE on the validation set 
## of each of the five folds.

### okay also muss man hier die kmeans() cluster funktion ändern, 
### und dann den Code von spatial cv nochmal ausführen
### der jetzt gerade noch nicht funktioniert.

dfs$clust_env <- as.factor(kmeans(
  dfs[,5:6],
  centers = 5
)$cluster)

# plot world map cluster distribution on mat+map
ggplot() +
  # plot coastline
  geom_sf(data = coast,
          colour = 'black',
          size = 0.2) +
  
  # set extent in longitude and latitude
  coord_sf(
    ylim = c(-60, 80),
    expand = FALSE) +  # to draw map strictly bounded by the specified extent
  
  # plot points on map
  geom_point(data = dfs, aes(x = lon, y = lat, color= clust_env), size = 0.6) +
  scale_color_manual(values=my_colors) +
  labs(x = "", y = "") +
  theme(legend.position = "bottom") +
  labs(title="Modelinput: mat+map")


## Plot the distribution of leaf N by cluster.
dfs %>%
  ggplot( aes(x=clust_env, y=leafN, fill=clust_env)) +
  geom_boxplot() +
  scale_fill_manual(values=my_colors) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")



dfs %>% 
  ggplot(aes(x=leafN)) +
  geom_histogram(fill = "brown3", color = "black") +
  facet_wrap(~clust_env) +
  labs(title = "Modelinput: mat+map")



## 2) Compare the results of the environmental cross-validation to the 
## results of the random and the spatial cross-validation and discuss reasons for 
## why you observe a difference in the cross-validation metrics (if you do).




