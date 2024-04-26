# load packages
suppressPackageStartupMessages({
  library("phyloseq")
  library("tidyverse")
  library("ggpubr")
  library("caret")
  library("knitr")
  library("mlbench")
  library("ranger")
  library("pheatmap")
  library("patchwork")
  library("ggplotify")
  library("pheatmap")
  
})

# set plotting theme
theme_set(theme_bw(14))
# load phyloseq object
ps1 <- readRDS("ps1.rds")

n_features <- nrow(otu_table(ps1))

##---------------------------------------------------------------------------
# model with entire sample and only ASV as predictors
# extract data from phyloseq object
set.seed(1390)
dataMatrix <-
  data.frame(PMI = sample_data(ps1)$PMI,
             snow_depth = sample_data(ps1)$snow_depth,
             t(otu_table(ps1)))

# create partition index
set.seed(1390)
trainIndex <- createDataPartition(dataMatrix$PMI,
                                  p = 0.7,
                                  list = FALSE,
                                  times = 1)

# data partition
set.seed(1390)
dataMatrixTrain <- dataMatrix[trainIndex,]
dataMatrixTest <- dataMatrix[-trainIndex,]

# drop unwanted viariables
drop <- "snow_depth"
dataMatrixTrain = dataMatrixTrain[,!(names(dataMatrixTrain) %in% drop)]
dataMatrixTest = dataMatrixTest[,!(names(dataMatrixTest) %in% drop)]

# fit initial ranger model 
set.seed(1390)
rf <- ranger(
  PMI ~ .,
  data = dataMatrixTrain,
  mtry = floor(n_features/3),
  respect.unordered.factors = "order"
)

# get OOB RMSE
(default_rmse <- sqrt(rf$prediction.error))

# # create hyperparameter grid
# hyper_grid <- expand.grid(
#   num.trees = floor(n_features / c(10, 20, 30, 40, 50)),
#   mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
#   min.node.size = c(1, 3, 5, 10),
#   replace = c(TRUE, FALSE),
#   sample.fraction = c(.5, .6, .7, .8, 1),
#   rmse = NA
# )

# set.seed(1390)
# #execute full cartesian grid search
# for (i in seq_len(nrow(hyper_grid))) {
#   # fit model for ith hyperparameter combination
#   fit <- ranger(
#     PMI ~ .,
#     data = dataMatrixTrain,
#     num.trees       = hyper_grid$num.trees[i],
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$min.node.size[i],
#     replace         = hyper_grid$replace[i],
#     sample.fraction = hyper_grid$sample.fraction[i],
#     verbose         = FALSE,
#     respect.unordered.factors = 'order',
#   )
#   # export OOB error
#   hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
# }
#
# # assess top models
# hyper_grid %>%
#   arrange(rmse) %>%
#   mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
#   head(5)

# num.trees mtry min.node.size replace sample.fraction     rmse perc_gain
# 1       417 1043             5   FALSE             0.8 2.779911  3.878964
# 2       417 1043             1   FALSE             0.8 2.786654  3.645805
# 3       417  626             1   FALSE             0.8 2.790947  3.497351
# 4       417 1389             3   FALSE             0.7 2.797370  3.275263
# 5       417 1389             1   FALSE             0.8 2.801655  3.127091

# build final model 
set.seed(1390)
rf_reg <- ranger(
  PMI ~ .,
  data = dataMatrixTrain,
  num.trees = 417,
  mtry = 1043,
  min.node.size = 1,
  sample.fraction = 0.8,
  importance = "permutation"
)

# model results
rf_reg

# error in train set
Metrics::rmse(dataMatrixTrain$PMI, rf_reg$predictions)
Metrics::mae(dataMatrixTrain$PMI, rf_reg$predictions)
# error in test set
Metrics::rmse(dataMatrixTest$PMI,
              predict(rf_reg, dataMatrixTest)$predictions)
Metrics::mae(dataMatrixTest$PMI,
             predict(rf_reg, dataMatrixTest)$predictions)

# Make predictions
predictions_train <- predict(rf_reg, dataMatrixTrain)$predictions
predictions_test <- predict(rf_reg, dataMatrixTest)$predictions
# create prediction dataframe
train_prediction <- data.frame(predictions_train, dataMatrixTrain$PMI)
test_prediction <- data.frame(predictions_test, dataMatrixTest$PMI)

# plot model
RF_tot <- ggplot() +
  geom_point(
    data = train_prediction,
    mapping = aes(x = dataMatrixTrain.PMI, y = predictions_train),
    color = '#440154',
    fill = '#443a83',
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = test_prediction,
    mapping = aes(x = dataMatrixTest.PMI, y = predictions_test),
    color = '#1fa187',
    fill = '#73d056',
    shape = 21,
    size = 3
  ) +
  theme_light(16) +
  geom_smooth(
    data = train_prediction,
    method = lm,
    mapping = aes(x = dataMatrixTrain.PMI, y = predictions_train),
    color = 'black',
    linetype = 2,
    se = FALSE
  ) +
  labs(title = "RF model - Entire sample", subtitle = "MAE = 2.47" ~ R ^
         2 ~ "= 0.78") +
  xlab("PMI (weeks)") + ylab("Estimated PMI (weeks)")

# normalised variable importance
imp_tot <- vip::vip(rf_reg, include_type = TRUE, scale=TRUE,
                    geom = "point", horizontal = TRUE, 
                    aesthetics = list(color = "#287D8EFF", shape = 16, size = 4)) +
  theme_light(16) +ggtitle("VIP") +
  labs(x = "Predictors", y = "Importance scores(%)",
       subtitle = "Total sample")  

# save plotted VIP
r <- rownames(tax_table(ps1)) %in% imp_tot$data$Variable
tax_table(ps1)[r,] %>%
  write.csv("vip.csv")

##---------------------------------------------------------------------------
# model including only swabs obtained from internal nose cavity and ASV as
# predictors
# subset internal sampling
ps1_int <-
  subset_samples(ps1, location == "Interior")

# select target variables
dataMatrix_int <-
  data.frame(
    PMI = sample_data(ps1_int)$PMI,
    snow_depth = sample_data(ps1_int)$snow_depth,
    t(otu_table(ps1_int))
  )

# set partition index
set.seed(1390)
trainIndex <- createDataPartition(dataMatrix_int$PMI,
                                  p = 0.7,
                                  list = FALSE,
                                  times = 1)

# perform splitting
set.seed(1390)
dataMatrixTrain_int <- dataMatrix_int[trainIndex, ]
dataMatrixTest_int <- dataMatrix_int[-trainIndex, ]

# drop variables not of interest
drop <- "snow_depth"
dataMatrixTrain_int = dataMatrixTrain_int[,!(names(dataMatrixTrain_int) %in% drop)]
dataMatrixTest_int = dataMatrixTest_int[,!(names(dataMatrixTest_int) %in% drop)]

# fit initial model 
set.seed(1390)
rf_int <- ranger(
  PMI ~ .,
  data = dataMatrixTrain_int,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

# get OOB RMSE
(default_rmse_int <- sqrt(rf_int$prediction.error))

# set.seed(1390)
# # execute full cartesian grid search
# for (i in seq_len(nrow(hyper_grid))) {
#   # fit model for ith hyperparameter combination
#   fit_int <- ranger(
#     PMI ~ .,
#     data = dataMatrixTrain_int,
#     num.trees       = hyper_grid$num.trees[i],
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$min.node.size[i],
#     replace         = hyper_grid$replace[i],
#     sample.fraction = hyper_grid$sample.fraction[i],
#     verbose         = FALSE,
#     respect.unordered.factors = 'order',
#   )
#   # export OOB error
#   hyper_grid$rmse[i] <- sqrt(fit_int$prediction.error)
# }
#
# # assess top models
# hyper_grid %>%
#   arrange(rmse) %>%
#   mutate(perc_gain = (default_rmse_int - rmse) / default_rmse * 100) %>%
#   head(5)
#
# num.trees mtry min.node.size replace sample.fraction
# 1       208 1043             5   FALSE             0.7
# 2        83 1669             3   FALSE             0.8
# 3       208 1043             3    TRUE             1.0
# 4       104 1669             3   FALSE             0.7
# 5       104 1389             3   FALSE             0.8
# rmse perc_gain
# 1 2.629805  2.875227
# 2 2.632166  2.793573
# 3 2.636006  2.660815
# 4 2.639843  2.528144
# 5 2.644397  2.370654

# fit final model
set.seed(1390)
rf_reg_int <- ranger(
  PMI ~ .,
  data = dataMatrixTrain_int,
  num.trees = 208,
  mtry = 625 ,
  min.node.size = 5,
  sample.fraction = .7,
  importance = "permutation"
)

# model results
rf_reg_int

# error in train set
Metrics::rmse(dataMatrixTrain_int$PMI, rf_reg_int$predictions)
Metrics::mae(dataMatrixTrain_int$PMI, rf_reg_int$predictions)
# error in test set
Metrics::rmse(dataMatrixTest_int$PMI,
              predict(rf_reg_int, dataMatrixTest_int)$predictions)
Metrics::mae(dataMatrixTest_int$PMI,
             predict(rf_reg_int, dataMatrixTest_int)$predictions)

# Make predictions
predictions_train_int <-
  predict(rf_reg_int, dataMatrixTrain_int)$predictions
predictions_test_int <-
  predict(rf_reg_int, dataMatrixTest_int)$predictions
# save prediction in dataframe
train_prediction_int <-
  data.frame(predictions_train_int, dataMatrixTrain_int$PMI)
test_prediction_int <-
  data.frame(predictions_test_int, dataMatrixTest_int$PMI)

# plot model
RF_int <- ggplot() +
  geom_point(
    data = train_prediction_int,
    mapping = aes(x = dataMatrixTrain_int.PMI, y = predictions_train_int),
    color = '#440154',
    fill = '#443a83',
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = test_prediction_int,
    mapping = aes(x = dataMatrixTest_int.PMI, y = predictions_test_int),
    color = '#1fa187',
    fill = '#73d056',
    shape = 21,
    size = 3
  ) +
  theme_light(16) +
  geom_smooth(
    data = train_prediction_int,
    method = lm,
    mapping = aes(x = dataMatrixTrain_int.PMI, y = predictions_train_int),
    color = 'black',
    linetype = 2,
    se = FALSE
  ) +
  labs(title = "RF model - Internal sample", subtitle = "MAE = 2.21," ~
         R ^ 2 ~ "= 0.81") +
  xlab("PMI (weeks)") + ylab("Estimated PMI (weeks)")

# normalised variable importance
imp_int <- vip::vip(rf_reg_int, include_type = TRUE, scale=TRUE,
                    geom = "point", horizontal = TRUE, 
                    aesthetics = list(color = "#287D8EFF", shape = 16, size = 4)) +
  theme_light(16) +ggtitle("VIP") +
  labs(x = "Predictors", y = "Importance scores (%)") +
  labs(x = "Predictors", y = "Importance scores(%)",
       subtitle = "Internal sample")  

# save plotted VIP
r <- rownames(tax_table(ps1_int)) %in% imp_int$data$Variable
tax_table(ps1_int)[r,] %>%
  write.csv("vip_int.csv")

##---------------------------------------------------------------------------
# model including only swabs obtained from external nose cavity and ASV as
# predictors
# subset external sample
ps1_ext <-
  subset_samples(ps1, location == "Exterior")

# select target variables
dataMatrix_ext <-
  data.frame(
    PMI = sample_data(ps1_ext)$PMI,
    snow_depth = sample_data(ps1_ext)$snow_depth,
    t(otu_table(ps1_ext))
  )

# set partition index
set.seed(1390)
trainIndex <- createDataPartition(dataMatrix_ext$PMI,
                                  p = 0.7,
                                  list = FALSE,
                                  times = 1)

# perform splitting
set.seed(1390)
dataMatrixTrain_ext <- dataMatrix_ext[trainIndex, ]
dataMatrixTest_ext <- dataMatrix_ext[-trainIndex, ]

# drop variables not of interest
drop <- "snow_depth"
dataMatrixTrain_ext = dataMatrixTrain_ext[,!(names(dataMatrixTrain_ext) %in% drop)]
dataMatrixTest_ext = dataMatrixTest_ext[,!(names(dataMatrixTest_ext) %in% drop)]

# fit initial model 
set.seed(1390)
rf_ext <- ranger(
  PMI ~ .,
  data = dataMatrixTrain_ext,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

# get OOB RMSE
(default_rmse_ext <- sqrt(rf_ext$prediction.error))

# set.seed(1930)
# # execute full cartesian grid search
# for (i in seq_len(nrow(hyper_grid))) {
#   # fit model for ith hyperparameter combination
#   fit_ext <- ranger(
#     PMI ~ .,
#     data = dataMatrixTrain_ext,
#     num.trees       = hyper_grid$num.trees[i],
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$min.node.size[i],
#     replace         = hyper_grid$replace[i],
#     sample.fraction = hyper_grid$sample.fraction[i],
#     verbose         = FALSE,
#     seed            = n_features/4,
#     respect.unordered.factors = 'order',
#   )
#   # export OOB error
#   hyper_grid$rmse[i] <- sqrt(fit_ext$prediction.error)
# }

# assess top models
# hyper_grid %>%
#   arrange(rmse) %>%
#   mutate(perc_gain = (default_rmse_ext - rmse) / default_rmse * 100) %>%
#   head(5)

# num.trees mtry min.node.size replace sample.fraction     rmse perc_gain
# 1        83  626             1   FALSE             0.8 3.141102  5.459538
# 2       104  626             1   FALSE             0.8 3.152205  5.075598
# 3       139  626             1   FALSE             0.8 3.174772  4.295332
# 4       208  626             1   FALSE             0.8 3.178509  4.166107
# 5       104 1389             1   FALSE             0.8 3.191688  3.710424

# fit final model
set.seed(1390)
rf_reg_ext <- ranger(
  PMI ~ .,
  data = dataMatrixTrain_ext,
  num.trees = 83,
  mtry =  626,
  min.node.size = 1,
  sample.fraction = .8,
  importance = "permutation"
)

# model results
rf_reg_ext

# Make predictions
predictions_train_ext <- predict(rf_reg_ext, dataMatrixTrain_ext)$predictions
predictions_test_ext <- predict(rf_reg_ext, dataMatrixTest_ext)$predictions
# save prediction in dataframe
train_prediction_ext <- data.frame(predictions_train_ext, dataMatrixTrain_ext$PMI)
test_prediction_ext <- data.frame(predictions_test_ext, dataMatrixTest_ext$PMI)

# error in train set
Metrics::rmse(dataMatrixTrain_ext$PMI, rf_reg_ext$predictions)
Metrics::mae(dataMatrixTrain_ext$PMI, rf_reg_ext$predictions)
# error in test set
Metrics::rmse(dataMatrixTest_ext$PMI, predict(rf_reg_ext, dataMatrixTest_ext)$predictions)
Metrics::mae(dataMatrixTest_ext$PMI, predict(rf_reg_ext, dataMatrixTest_ext)$predictions)

# plot model
RF_ext <- ggplot() +
  geom_point(
    data = train_prediction_ext,
    mapping = aes(x = dataMatrixTrain_ext.PMI, y = predictions_train_ext),
    color = '#440154',
    fill = '#443a83',
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = test_prediction_ext,
    mapping = aes(x = dataMatrixTest_ext.PMI, y = predictions_test_ext),
    color = '#1fa187',
    fill = '#73d056',
    shape = 21,
    size = 3
  ) +
  theme_light(16) +
  geom_smooth(
    data = train_prediction_ext,
    method = lm,
    mapping = aes(x = dataMatrixTrain_ext.PMI, y = predictions_train_ext),
    color = 'black',
    linetype = 2,
    se = FALSE
  ) +
  labs(title = "RF model - External sample", subtitle = "MAE = 2.89," ~
         R ^ 2 ~ "= 0.73") +
  xlab("PMI (weeks)") + ylab("Estimated PMI (weeks)")

# normalised variable importance
imp_ext <- vip::vip(rf_reg_ext, include_type = TRUE, scale=TRUE,
                    geom = "point", horizontal = TRUE, 
                    aesthetics = list(color = "#287D8EFF", shape = 16, size = 4)) +
  theme_light(16) +ggtitle("VIP") +
  labs(x = "Predictors", y = "Importance scores (%)") +
  labs(x = "Predictors", y = "Importance scores(%)",
       subtitle = "External sample") 

# save plotted VIP
r <- rownames(tax_table(ps1_ext)) %in% imp_ext$data$Variable
tax_table(ps1_ext)[r,] %>% 
  write.csv("vip_ext.csv")

#-------------------------------------------------------------------------------
# model including entire sample with temperature and snow coverage
set.seed(1390)

# select target variables
df <-
  data.frame(PMI = sample_data(ps1)$PMI,
             snow_depth = sample_data(ps1)$snow_depth,
             Temp = sample_data(ps1)$T,
             t(otu_table(ps1)))

# set partition index
set.seed(1390)
trainIndex <- createDataPartition(df$PMI,
                                  p = 0.7,
                                  list = FALSE,
                                  times = 1)

# perform splitting
set.seed(1390)
dfTrain <- df[trainIndex, ]
dfTest <- df[-trainIndex, ]

# fit initial model 
set.seed(1390)
rf <- ranger(
  PMI ~ .,
  data = dfTrain,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

# get OOB RMSE
(default_rmse <- sqrt(rf$prediction.error))

# # create hyperparameter grid
# hyper_grid <- expand.grid(
#   num.trees = floor(n_features / c(10, 20, 30, 40, 50)),
#   mtry = floor(n_features * c(.05, .15, .25, .333, .4)),
#   min.node.size = c(1, 3, 5, 10),
#   replace = c(TRUE, FALSE),
#   sample.fraction = c(.5, .6, .7, .8, 1),
#   rmse = NA
# )
# 
# set.seed(1390)
# #execute full cartesian grid search
# for (i in seq_len(nrow(hyper_grid))) {
#   # fit model for ith hyperparameter combination
#   fit <- ranger(
#     PMI ~ .,
#     data = dfTrain,
#     num.trees       = hyper_grid$num.trees[i],
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$min.node.size[i],
#     replace         = hyper_grid$replace[i],
#     sample.fraction = hyper_grid$sample.fraction[i],
#     verbose         = FALSE,
#     respect.unordered.factors = 'order',
#   )
#   # export OOB error
#   hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
# }
# 
# # assess top models
# hyper_grid %>%
#   arrange(rmse) %>%
#   mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
#   head(5)
# 
# num.trees mtry min.node.size replace sample.fraction     rmse perc_gain
# 1       139 1669             1   FALSE             0.8 1.606841 12.730570
# 2       208 1669             1   FALSE             0.8 1.612837 12.404922
# 3       417 1669             5   FALSE             0.8 1.647693 10.511812
# 4       139 1669             3   FALSE             0.8 1.666822  9.472911
# 5       417 1669             1   FALSE             0.8 1.670299  9.28409

# fit final model
set.seed(1390)
rf_reg <- ranger(
  PMI ~ .,
  data = dfTrain,
  num.trees = 139,
  mtry = 1669,
  min.node.size = 1,
  sample.fraction = 0.8,
  importance = "permutation"
)

# model results
rf_reg

# error in train set
Metrics::rmse(dfTrain$PMI, rf_reg$predictions)
Metrics::mae(dfTrain$PMI, rf_reg$predictions)
# error in test set
Metrics::rmse(dfTest$PMI, predict(rf_reg, dfTest)$predictions)
Metrics::mae(dfTest$PMI, predict(rf_reg, dfTest)$predictions)

# Make predictions
predictions_train <- predict(rf_reg, dfTrain)$predictions
predictions_test <- predict(rf_reg, dfTest)$predictions
# save prediction in dataframe
train_prediction <- data.frame(predictions_train, dfTrain$PMI)
test_prediction <- data.frame(predictions_test, dfTest$PMI)

#  plot model
RF_tot_complex <- ggplot() +
  geom_point(
    data = train_prediction,
    mapping = aes(x = dfTrain.PMI, y = predictions_train),
    color = '#440154',
    fill = '#443a83',
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = test_prediction,
    mapping = aes(x = dfTest.PMI, y = predictions_test),
    color = '#1fa187',
    fill = '#73d056',
    shape = 21,
    size = 3
  ) +
  theme_light(16) +
  geom_smooth(
    data = train_prediction,
    method = lm,
    mapping = aes(x = dfTrain.PMI, y = predictions_train),
    color = 'black',
    linetype = 2,
    se = FALSE
  ) +
  labs(title = "RF model - Entire sample", subtitle = "MAE = 1.36," ~R^2~"= 0.91") +
  xlab("PMI (weeks)") + ylab("Estimated PMI (weeks)")

# normalised variable importance
imp_complex <- vip::vip(rf_reg, include_type = TRUE, scale=TRUE,
                            geom = "point", horizontal = TRUE, 
                            aesthetics = list(color = "#443a83", shape = 16, size = 4)) +
  theme_light(16) +ggtitle("VIP") +
  labs(x = "Predictors", y = "Importance scores(%)",
       subtitle = "Total complex") 

# save plotted VIP
r <- rownames(tax_table(ps1_int)) %in% imp_complex$data$Variable
tax_table(ps1_int)[r,] %>% 
  write.csv("vip_complex.csv")

##---------------------------------------------------------------------------
# model including only swabs obtained from internal nose cavity, ASV, 
# temperature, and snow coverage as predictors
# subset internal sampling
ps1_int <-
  subset_samples(ps1, location == "Interior")

# select target variables
df_int <-
  data.frame(
    PMI = sample_data(ps1_int)$PMI,
    snow_depth = sample_data(ps1_int)$snow_depth,
    Temp = sample_data(ps1_int)$T,
    t(otu_table(ps1_int))
  )

# set partition index
set.seed(1390)
trainIndex <- createDataPartition(df_int$PMI,
                                  p = 0.7,
                                  list = FALSE,
                                  times = 1)

# perform splitting
set.seed(1390)
dfTrain_int <- df_int[trainIndex,]
dfTest_int <- df_int[-trainIndex,]

# fit initial model 
set.seed(1390)
rf_int <- ranger(
  PMI ~ .,
  data = dfTrain_int,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

# get OOB RMSE
(default_rmse_int <- sqrt(rf_int$prediction.error))

# # execute full cartesian grid search
# for (i in seq_len(nrow(hyper_grid))) {
#   # fit model for ith hyperparameter combination
#   fit_int <- ranger(
#     PMI ~ .,
#     data = dfTrain_int,
#     num.trees       = hyper_grid$num.trees[i],
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$min.node.size[i],
#     replace         = hyper_grid$replace[i],
#     sample.fraction = hyper_grid$sample.fraction[i],
#     verbose         = FALSE,
#     respect.unordered.factors = 'order',
#   )
#   # export OOB error
#   hyper_grid$rmse[i] <- sqrt(fit_int$prediction.error)
# }
# 
# # assess top models
# hyper_grid %>%
#   arrange(rmse) %>%
#   mutate(perc_gain = (default_rmse_int - rmse) / default_rmse * 100) %>%
#   head(5)


# num.trees mtry min.node.size replace sample.fraction     rmse perc_gain
# 1        83 1389             5   FALSE             0.8 1.949741  8.295515
# 2       139 1669             5   FALSE             0.7 1.988829  6.172612
# 3        83 1389             1   FALSE             0.8 2.009306  5.060448
# 4       208 1389             3   FALSE             0.8 2.022049  4.368364
# 5       208 1669             3   FALSE             0.8 2.024367  4.242486

# fit initial model 
set.seed(1390)
rf_reg_int <- ranger(
  PMI ~ .,
  data = dfTrain_int,
  num.trees = 83,
  mtry = 1389 ,
  min.node.size = 5,
  sample.fraction = .8,
  importance = "permutation"
)

# model results
rf_reg_int

# error in train set
Metrics::rmse(dfTrain_int$PMI, rf_reg_int$predictions)
Metrics::mae(dfTrain_int$PMI, rf_reg_int$predictions)
# error in test set
Metrics::rmse(dfTest_int$PMI, predict(rf_reg_int, dfTest_int)$predictions)
Metrics::mae(dfTest_int$PMI, predict(rf_reg_int, dfTest_int)$predictions)

# Make predictions
predictions_train_int <- predict(rf_reg_int, dfTrain_int)$predictions
predictions_test_int <- predict(rf_reg_int, dfTest_int)$predictions
# save prediction in dataframe
train_prediction_int <- data.frame(predictions_train_int, dfTrain_int$PMI)
test_prediction_int <- data.frame(predictions_test_int, dfTest_int$PMI)

# plot model
RF_int_complex <- ggplot() +
  geom_point(
    data = train_prediction_int,
    mapping = aes(x = dfTrain_int.PMI, y = predictions_train_int),
    color = '#440154',
    fill = '#443a83',
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = test_prediction_int,
    mapping = aes(x = dfTest_int.PMI, y = predictions_test_int),
    color = '#1fa187',
    fill = '#73d056',
    shape = 21,
    size = 3
  ) +
  theme_light(16) +
  geom_smooth(
    data = train_prediction_int,
    method = lm,
    mapping = aes(x = dfTrain_int.PMI, y = predictions_train_int),
    color = 'black',
    linetype = 2,
    se = FALSE
  ) +
  labs(title = "RF model - Internal sample", subtitle = "MAE = 1.60," ~R^2~"= 0.88") +
  xlab("PMI (weeks)") + ylab("Estimated PMI (weeks)")

# normalised variable importance
imp_int_complex <- vip::vip(rf_reg_int, include_type = TRUE, scale=TRUE,
                        geom = "point", horizontal = TRUE, 
                        aesthetics = list(color = "#443a83", shape = 16, size = 4)) +
  theme_light(16) +ggtitle("VIP") +
  labs(x = "Predictors", y = "Importance scores(%)",
       subtitle = "Internal complex") 

# save plotted VIP
r <- rownames(tax_table(ps1_int)) %in% imp_int_complex$data$Variable
tax_table(ps1_int)[r,] %>% 
  write.csv("vip_int_complex.csv")

##---------------------------------------------------------------------------
# model including only swabs obtained from external nose cavity, ASV, 
# temperature, and snow coverage as predictors
# subset external sampling
ps1_ext <-
  subset_samples(ps1, location == "Exterior")

# select target variables
df_ext <-
  data.frame(
    PMI = sample_data(ps1_ext)$PMI,
    snow_depth = sample_data(ps1_ext)$snow_depth,
    Temp = sample_data(ps1_ext)$T,
    t(otu_table(ps1_ext))
  )

# set partition index
set.seed(1390)
trainIndex <- createDataPartition(df_ext$PMI,
                                  p = 0.7,
                                  list = FALSE,
                                  times = 1)


# perform splitting
set.seed(1390)
dfTrain_ext <- df_ext[trainIndex,]
dfTest_ext <- df_ext[-trainIndex,]

# fit initial model 
set.seed(1390)
rf_ext <- ranger(
  PMI ~ .,
  data = dfTrain_ext,
  mtry = floor(n_features / 3),
  respect.unordered.factors = "order"
)

# get OOB RMSE
(default_rmse_ext <- sqrt(rf_ext$prediction.error))

# set.seed(n_features/3)
# # execute full cartesian grid search
# for (i in seq_len(nrow(hyper_grid))) {
#   # fit model for ith hyperparameter combination
#   fit_ext <- ranger(
#     PMI ~ .,
#     data = dfTrain_ext,
#     num.trees       = hyper_grid$num.trees[i],
#     mtry            = hyper_grid$mtry[i],
#     min.node.size   = hyper_grid$min.node.size[i],
#     replace         = hyper_grid$replace[i],
#     sample.fraction = hyper_grid$sample.fraction[i],
#     verbose         = FALSE,
#     seed            = n_features/4,
#     respect.unordered.factors = 'order',
#   )
#   # export OOB error
#   hyper_grid$rmse[i] <- sqrt(fit_ext$prediction.error)
# }
# 
# # assess top models
# hyper_grid %>%
#   arrange(rmse) %>%
#   mutate(perc_gain = (default_rmse_ext - rmse) / default_rmse * 100) %>%
#   head(5)

# num.trees mtry min.node.size replace sample.fraction     rmse perc_gain
# 1       104 1389             3   FALSE             0.8 2.202808  15.26774
# 2       417 1389             1   FALSE             0.8 2.214940  14.60883
# 3        83 1389             3   FALSE             0.8 2.230122  13.78427
# 4       104 1389             5   FALSE             0.8 2.231860  13.68990
# 5       208 1389             1   FALSE             0.8 2.235123  13.51266

# fit final model
set.seed(1390)
rf_reg_ext <- ranger(
  PMI ~ .,
  data = dfTrain_ext,
  num.trees = 104,
  mtry =  1389,
  min.node.size = 3,
  sample.fraction = .8,
  importance = "permutation"
)

# model results
rf_reg_ext

# Make predictions
predictions_train_ext <- predict(rf_reg_ext, dfTrain_ext)$predictions
predictions_test_ext <- predict(rf_reg_ext, dfTest_ext)$predictions
# save prediction in dataframe
train_prediction_ext <- data.frame(predictions_train_ext, dfTrain_ext$PMI)
test_prediction_ext <- data.frame(predictions_test_ext, dfTest_ext$PMI)

# error in train set
Metrics::rmse(dfTrain_ext$PMI, rf_reg_ext$predictions)
Metrics::mae(dfTrain_ext$PMI, rf_reg_ext$predictions)
# error in test set
Metrics::rmse(dfTest_ext$PMI, predict(rf_reg_ext, dfTest_ext)$predictions)
Metrics::mae(dfTest_ext$PMI, predict(rf_reg_ext, dfTest_ext)$predictions)

# plot model
RF_ext_complex <- ggplot() +
  geom_point(
    data = train_prediction_ext,
    mapping = aes(x = dfTrain_ext.PMI, y = predictions_train_ext),
    color = '#440154',
    fill = '#443a83',
    shape = 24,
    size = 3
  ) +
  geom_point(
    data = test_prediction_ext,
    mapping = aes(x = dfTest_ext.PMI, y = predictions_test_ext),
    color = '#1fa187',
    fill = '#73d056',
    shape = 21,
    size = 3
  ) +
  theme_light(16) +
  geom_smooth(
    data = train_prediction_ext,
    method = lm,
    mapping = aes(x = dfTrain_ext.PMI, y = predictions_train_ext),
    color = 'black',
    linetype = 2,
    se = FALSE
  ) +
  labs(title = "RF model - External sample", subtitle = "MAE = 1.83," ~R^2~"= 0.87") +
  xlab("PMI (weeks)") + ylab("Estimated PMI (weeks)")

# normalised variable importance
imp_ext_complex <- vip::vip(rf_reg_ext, include_type = TRUE, scale=TRUE,
                            geom = "point", horizontal = TRUE, 
                            aesthetics = list(color = "#443a83", shape = 16, size = 4)) +
  theme_light(16) +ggtitle("VIP") +
  labs(x = "Predictors", y = "Importance scores(%)",
       subtitle = "Internal complex") 

# save plotted VIP
r <- rownames(tax_table(ps1_int)) %in% imp_ext_complex$data$Variable
tax_table(ps1_int)[r,] %>% 
  write.csv("vip_ext_complex.csv")

#-------------------------------------------------------------------------------
# scatterplot of all predictions
figure_complex <- ggarrange(
  RF_tot,
  RF_int,
  RF_ext,
  RF_tot_complex,
  RF_int_complex,
  RF_ext_complex,
  ncol = 3,
  nrow = 2,
  common.legend = TRUE,
  labels = 'auto'
)

# annotate figre
annotate_figure(
  figure_complex,
  bottom = text_grob(
    "Triangle: Train set; Circle: Test set",
    color = "black",
    hjust = 1.05,
    x = 1,
    face = "bold",
    size = 14
  )
  
)

# save plot
ggsave("figure6.pdf", width = 14.98, height = 10.36)

##---------------------------------------------------------------------------
# prepare heatmap for VIP across PMI 
# prepare data
imps1 <- data.frame(var = names(dfTrain)[-1],
                    imps1 = rf_reg$variable.importance)
imp.sort <- arrange(imps1, desc(imps1))
imp.25 <- imp.sort[imp.sort$imps1 > 1, ]

r <- rownames(tax_table(ps1)) %in% imp.25$var
tax_table(ps1)[r, ]

heat_subset <-
  subset(otu_table(ps1), rownames(otu_table(ps1)) %in% c(rownames(otu_table(ps1)[r,])))
new_physeq <-
  merge_phyloseq(heat_subset, tax_table(ps1), sample_data(ps1))

heat <-
  data.frame(PMI = sample_data(new_physeq)$PMI,
             t(otu_table(new_physeq)))

group_mean <- aggregate(x = heat,
                        by = list(heat$PMI),
                        FUN = mean)

group_mean <- group_mean[-1]
rownames(group_mean) <- group_mean$PMI
group_mean <- group_mean[-1]

heatmap <- 
  pheatmap::pheatmap(
    scale(t(group_mean)),
    cluster_cols = FALSE,
    color = viridis::viridis(50),
    angle_col = 45,
    fontsize_row = 14,
    fontsize_col = 14,
    main = "RF Complex VIP - Total Sample"
  )

# combine plots
patchwork <- (imp_tot | imp_int | imp_ext) /
  (imp_complex | imp_int_complex | imp_ext_complex) /
  as.ggplot(heatmap) + 
  plot_layout(nrow = 3, widths = c(0.5, 0.5, 0.5), 
              heights = unit(c(4, 4, 9), c('cm', 'cm', 'cm')))

# combine plots
patchwork +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 14))

# save plot 
ggsave("figure7.pdf", width = 11.34, height = 11.28)


