
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
  
  rmse <- sqrt(mean(lm(pred$predictions~df$leafN[idx_val)$residuals^2)) 
  # the root mean square error on the validation set
  
  return(tibble(rsq = rsq, rmse = rmse))
}

# apply function on each custom fold and collect validation results in a nice
# data frame
out <- purrr::map2_dfr(
  group_folds_train,
  group_folds_test,
  ~train_test_by_fold(dfs, group_folds_train, group_folds_test) # das muss man hier definieren
) |> 
  mutate(test_fold = 1:5)
out
