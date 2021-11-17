### Supplemental Materials - 1 (Data and scripts)
### Main modelling and anomaly detection

## ---- load
library(tidyverse)
library(lubridate)
library(ggpubr)
library(forecast)
library(visreg)
library(gridExtra)
library(conduits)
library(mgcv)
library(forecast)
library(imputeTS)

## loadExternalFunctions
source(here::here("presentation",
                  "scripts", "SPOT_algorithm.R"))

## loadPRINData
load(here::here("code-sharing", "data", "data_lagged_upstream.rda"))
load(here::here("code-sharing", "data", "data_downstream_imputed.rda"))

load(here::here("code-sharing", "data", "PRIN_5min_flagged.rda"))
# load(here::here("weird_river_data_paper", "data", "PRIN_5min_cleaned.rda"))



## ---- dataPreparation
# taking flagged data for turbidity
turbidity_downstream_flagged <- PRIN_5min_flagged %>%
  rename(Timestamp = roundedTimestamp) %>%
  mutate(Timestamp = ymd_hms(Timestamp)) %>% 
  filter(site == "down", Timestamp >= ymd("2019-10-01") 
         & Timestamp < ymd("2020-01-01")) %>% 
  select(Timestamp, turbidity, turbidityAnomalyFlag) %>% 
  rename("turbidity_downstream" = turbidity)
    
# other downstream variables from the cleaned imputed data
data_downstream <- data_downstream_imputed %>% 
  select(Timestamp, level_downstream, temperature_downstream,
         conductance_downstream, dissolved_oxygen_downstream,
         pH_downstream)
  
  
# combining all data
data_all <- list(data_downstream, data_lagged_upstream, 
                 turbidity_downstream_flagged) %>% 
  purrr::reduce(left_join, by = "Timestamp")


## scaling-and-transformations

# computing log of turbidity
data_all <- data_all %>%
  mutate(turbidity_upstream_dt_log = log(turbidity_upstream_dt),
         turbidity_downstream_log = log(turbidity_downstream),
         turbidity_downstream_log_lag1 = lag(turbidity_downstream_log, 1),
         turbidity_downstream_log_lag2 = lag(turbidity_downstream_log, 2),
         turbidity_downstream_log_lag3 = lag(turbidity_downstream_log, 3))


## training-data
# taking data from "2019-10-01" to "2019-11-30" as training data
data_train <- data_all %>%
  filter(Timestamp >= ymd("2019-10-01")
         & Timestamp < ymd("2019-12-01"))


# We will remove the anomalies in the turbidity_downstream in training data
data_train <- data_train %>%
  mutate(turbidity_downstream = if_else(turbidityAnomalyFlag=="outlier",
                                        as.numeric(NA),
                                        turbidity_downstream))


## ---- Model-fitting

######################
# Fitting-GAM-up-AR
######################

#taking level, temperature, turbidity and conductance from the upstream
#sensor as predictors in this model
gam_up_AR <- gam(turbidity_downstream_log ~
                   s(turbidity_downstream_log_lag1) +
                   s(turbidity_downstream_log_lag2) +
                   s(turbidity_downstream_log_lag3) +
                   s(turbidity_upstream_dt_log) +
                   s(level_upstream_dt) +
                   s(temperature_upstream_dt) +
                   s(conductance_upstream_dt),
                  data = data_train)


## compute-resids-GAM-up-AR

data_residuals <- data_all %>%
  select(Timestamp, turbidity_downstream, turbidity_downstream_log,
         turbidityAnomalyFlag) %>%
  mutate(predict_mod_up_AR = as.numeric(predict.gam(gam_up_AR, newdata = data_all)),
         residuals_up_AR = turbidity_downstream_log - predict_mod_up_AR,
         turbidityAnomalyFlag = case_when(turbidityAnomalyFlag == 0 ~ "typical",
                                          turbidityAnomalyFlag == 1 ~ "outlier",
                                          TRUE ~ as.character(NA)))



## SPOT-for-GAM-up-AR

#detecting outliers in residuals from the model using SPOT algorithm

q <- 0.01
t_prob <- 0.98
n = floor(nrow(data_residuals)*0.30) # first 30% of the data (might contain missing values)


# Using SPOT and computing confusion matrix
confusion_mod_up_AR <- compute_confusion(data = data_residuals,
                                         residuals = "residuals_up_AR",
                                         name = "model_up_AR", n = n,
                                         q = q, t_prob = t_prob,
                                         cmax = FALSE)

confusion_mod_up_AR <- confusion_mod_up_AR %>%
  rename(AnomalyFlag_SPOT_up_AR = AnomalyFlag_SPOT,
         confusion_mod_up_AR = confusion_matrix)

data_residuals <- data_residuals %>%
  left_join(confusion_mod_up_AR, by = "Timestamp")


## classification-GAM-up-AR

classif_plot_turb_up_AR <- data_residuals %>%
  select(Timestamp, turbidity_downstream, confusion_mod_up_AR) %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(Timestamp, turbidity_downstream,
                 color = confusion_mod_up_AR,
                 shape = confusion_mod_up_AR,
                 size = confusion_mod_up_AR,
                 alpha = confusion_mod_up_AR)) +
  scale_color_manual(values = c("#009E73", "#999999", "#E69F00",
                                "#FF3333")) +
  scale_shape_manual(values = c(17,16,18,15)) +
  scale_size_manual(values = c(2,0.3,2,2)) +
  scale_alpha_manual(values = c(1,1,1,1)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ylab("Turbidity (FNU)")

calib_time <- data_residuals %>% 
  slice(1:n) %>% 
  pull(Timestamp)

TP_time <- data_residuals %>% 
  filter(confusion_mod_up_AR == "TP") %>% 
  pull(Timestamp)

TN_time <- data_residuals %>% 
  filter(confusion_mod_up_AR == "TN") %>% 
  pull(Timestamp)

FP_time <- data_residuals %>% 
  filter(confusion_mod_up_AR == "FP") %>% 
  pull(Timestamp)

FN_time <- data_residuals %>% 
  filter(confusion_mod_up_AR == "FN") %>% 
  pull(Timestamp)

classif_plot_resid_up_AR <- data_residuals %>%
  select(Timestamp, residuals_up_AR) %>%
  drop_na() %>%
  ggplot() +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% calib_time),
             aes(Timestamp, residuals_up_AR), 
             color = "#000000",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp > max(calib_time) &
                        Timestamp %in% TN_time),
             aes(Timestamp, residuals_up_AR),
             color = "#999999",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% TP_time),
             aes(Timestamp, residuals_up_AR),
             color = "#009E73",
             shape = 17,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FP_time),
             aes(Timestamp, residuals_up_AR),
             color = "#E69F00",
             shape = 18,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FN_time),
             aes(Timestamp, residuals_up_AR),
             color = "#FF3333",
             shape = 15,
             size = 2) +
  ylab("residuals")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(classif_plot_turb_up_AR)


######################
  ## Fitting-GAM-up
######################

gam_up <- gam(turbidity_downstream_log ~
                s(turbidity_upstream_dt_log) +
                s(dissolved_oxygen_upstream_dt) +
                s(conductance_upstream_dt) +
                s(level_upstream_dt) +
                s(temperature_upstream_dt),
              data = data_train)


## SPOT-for-GAM-up

# computing residuals using the Model-up from all data

data_residuals <- data_residuals %>%
  mutate(predict_mod_up = as.numeric(predict.gam(gam_up, newdata = data_all)),
         residuals_up= turbidity_downstream_log - predict_mod_up)

#choosing the same hyper-parameters as GAM-up-AR residuals
confusion_mod_up <- compute_confusion(data = data_residuals,
                                         residuals = "residuals_up",
                                         name = "model_up", n = n,
                                         q = q, t_prob = t_prob,
                                         cmax = FALSE)

confusion_mod_up <- confusion_mod_up %>%
  rename(AnomalyFlag_SPOT_up = AnomalyFlag_SPOT,
         confusion_mod_up = confusion_matrix)
data_residuals <- data_residuals %>%
  left_join(confusion_mod_up, by = "Timestamp")



## classification-GAM-up

classif_plot_turb_up <- data_residuals %>%
  select(Timestamp, turbidity_downstream, confusion_mod_up) %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(Timestamp, turbidity_downstream,
                 color = confusion_mod_up,
                 shape = confusion_mod_up,
                 size = confusion_mod_up,
                 alpha = confusion_mod_up)) +
  scale_color_manual(values = c("#009E73", "#999999", "#E69F00",
                                "#FF3333")) +
  scale_shape_manual(values = c(17,16,18,15)) +
  scale_size_manual(values = c(2,0.3,2,2)) +
  scale_alpha_manual(values = c(1,1,1,1)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ylab("Turbidity (FNU)")

calib_time <- data_residuals %>% 
  slice(1:n) %>% 
  pull(Timestamp)

TP_time <- data_residuals %>% 
  filter(confusion_mod_up == "TP") %>% 
  pull(Timestamp)

TN_time <- data_residuals %>% 
  filter(confusion_mod_up == "TN") %>% 
  pull(Timestamp)

FP_time <- data_residuals %>% 
  filter(confusion_mod_up == "FP") %>% 
  pull(Timestamp)

FN_time <- data_residuals %>% 
  filter(confusion_mod_up == "FN") %>% 
  pull(Timestamp)

classif_plot_resid_up <- data_residuals %>%
  select(Timestamp, residuals_up) %>%
  drop_na() %>%
  ggplot() +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% calib_time),
             aes(Timestamp, residuals_up), 
             color = "#000000",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp > max(calib_time) &
                        Timestamp %in% TN_time),
             aes(Timestamp, residuals_up),
             color = "#999999",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% TP_time),
             aes(Timestamp, residuals_up),
             color = "#009E73",
             shape = 17,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FP_time),
             aes(Timestamp, residuals_up),
             color = "#E69F00",
             shape = 18,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FN_time),
             aes(Timestamp, residuals_up),
             color = "#FF3333",
             shape = 15,
             size = 2) +
  ylab("residuals")


########################
  ## Fitting-GAM-down
########################

#taking level, temperature, conductance and dissolved oxygen at lag1 from the
#downstream sensor as predictors in this model

# #computing lag1 of predictors
# data_train <- data_train %>%
#   mutate(level_downstream_lag1 = lag(level_downstream),
#          temperature_downstream_lag1 = lag(temperature_downstream),
#          conductance_downstream_lag1 = lag(conductance_downstream),
#          dissolved_oxygen_downstream_lag1 = lag(dissolved_oxygen_downstream))


gam_down <- gam(turbidity_downstream_log ~
                  s(conductance_downstream) +
                  s(dissolved_oxygen_downstream) +
                  s(level_downstream, k=6) +
                  s(temperature_downstream),
              data = data_train)



## SPOT-for-GAM-down

# computing residuals using the Model_B from all data

# data_all <- data_all %>%
#   mutate(level_downstream_lag1 = lag(level_downstream),
#          temperature_downstream_lag1 = lag(temperature_downstream),
#          conductance_downstream_lag1 = lag(conductance_downstream),
#          dissolved_oxygen_downstream_lag1 = lag(dissolved_oxygen_downstream))

data_residuals <- data_residuals %>%
  mutate(predict_mod_down = as.numeric(predict.gam(gam_down, newdata = data_all)),
         residuals_down = turbidity_downstream_log - predict_mod_down)

confusion_mod_down <- compute_confusion(data = data_residuals,
                                      residuals = "residuals_down",
                                      name = "model_down", n = n,
                                      q = q, t_prob = t_prob,
                                      cmax = FALSE)

confusion_mod_down <- confusion_mod_down %>%
  rename(AnomalyFlag_SPOT_down = AnomalyFlag_SPOT,
         confusion_mod_down = confusion_matrix)
data_residuals <- data_residuals %>%
  left_join(confusion_mod_down, by = "Timestamp")


## classification-GAM-down

classif_plot_turb_down <- data_residuals %>%
  select(Timestamp, turbidity_downstream, confusion_mod_down) %>%
  mutate(confusion_mod_down = factor(confusion_mod_down,
                            levels = c("TP", "TN", "FP", "FN",
                                       "NA"))) %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(Timestamp, turbidity_downstream,
                 color = confusion_mod_down,
                 shape = confusion_mod_down,
                 size = confusion_mod_down,
                 alpha = confusion_mod_down)) +
  scale_color_manual(values = c("#009E73", "#999999",
                                "#FF3333")) +
  scale_shape_manual(values = c(17,16,18,15)) +
  scale_size_manual(values = c(2,0.3,2,2)) +
  scale_alpha_manual(values = c(1,1,1,1)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ylab("Turbidity (FNU)")

calib_time <- data_residuals %>% 
  slice(1:n) %>% 
  pull(Timestamp)

TP_time <- data_residuals %>% 
  filter(confusion_mod_down == "TP") %>% 
  pull(Timestamp)

TN_time <- data_residuals %>% 
  filter(confusion_mod_down == "TN") %>% 
  pull(Timestamp)

FP_time <- data_residuals %>% 
  filter(confusion_mod_down == "FP") %>% 
  pull(Timestamp)

FN_time <- data_residuals %>% 
  filter(confusion_mod_down == "FN") %>% 
  pull(Timestamp)


classif_plot_resid_down <- data_residuals %>%
  select(Timestamp, residuals_down) %>%
  drop_na() %>%
  ggplot() +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% calib_time),
             aes(Timestamp, residuals_down), 
             color = "#000000",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp > max(calib_time) &
                        Timestamp %in% TN_time),
             aes(Timestamp, residuals_down),
             color = "#999999",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% TP_time),
             aes(Timestamp, residuals_down),
             color = "#009E73",
             shape = 17,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FP_time),
             aes(Timestamp, residuals_down),
             color = "#E69F00",
             shape = 18,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FN_time),
             aes(Timestamp, residuals_down),
             color = "#FF3333",
             shape = 15,
             size = 2) +
  ylab("residuals")


###########################
  ## Fitting-GAM-up-down
###########################

#adding temperature and conductance at lag1 from the downstream sensor 
#to the GAM-up

gam_up_down <- update(gam_up, .~.+
                   s(conductance_downstream) +
                   s(temperature_downstream))



## SPOT-for-GAM-up-down

# computing residuals using the model 3 from all data

data_residuals <- data_residuals %>%
  mutate(predict_mod_up_down = as.numeric(predict.gam(gam_up_down, newdata = data_all)),
         residuals_up_down = turbidity_downstream_log - predict_mod_up_down)


confusion_mod_up_down <- compute_confusion(data = data_residuals, 
                                        residuals = "residuals_up_down",
                                        name = "model_up_down", n = n,
                                        q = q, t_prob = t_prob,
                                        cmax = FALSE)

confusion_mod_up_down <- confusion_mod_up_down %>%
  rename(AnomalyFlag_SPOT_up_down = AnomalyFlag_SPOT,
         confusion_mod_up_down = confusion_matrix)
data_residuals <- data_residuals %>%
  left_join(confusion_mod_up_down, by = "Timestamp")



## classification-GAM-up-down

classif_plot_turb_up_down <- data_residuals %>%
  select(Timestamp, turbidity_downstream, confusion_mod_up_down) %>%
  mutate(confusion_mod_up_down = factor(confusion_mod_up_down,
                                        levels = c("TP", "TN", "FP", "FN",
                                                   "NA"))) %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(Timestamp, turbidity_downstream,
                 color = confusion_mod_up_down,
                 shape = confusion_mod_up_down,
                 size = confusion_mod_up_down,
                 alpha = confusion_mod_up_down)) +
  scale_color_manual(values = c("#009E73", "#999999", "#E69F00",
                                "#FF3333")) +
  scale_shape_manual(values = c(17,16,18,15)) +
  scale_size_manual(values = c(2,0.3,2,2)) +
  scale_alpha_manual(values = c(1,1,1,1)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ylab("Turbidity (FNU)")

calib_time <- data_residuals %>% 
  slice(1:n) %>% 
  pull(Timestamp)

TP_time <- data_residuals %>% 
  filter(confusion_mod_up_down == "TP") %>% 
  pull(Timestamp)

TN_time <- data_residuals %>% 
  filter(confusion_mod_up_down == "TN") %>% 
  pull(Timestamp)

FP_time <- data_residuals %>% 
  filter(confusion_mod_up_down == "FP") %>% 
  pull(Timestamp)

FN_time <- data_residuals %>% 
  filter(confusion_mod_up_down == "FN") %>% 
  pull(Timestamp)


classif_plot_resid_up_down <- data_residuals %>%
  select(Timestamp, residuals_up_down) %>%
  drop_na() %>%
  ggplot() +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% calib_time),
             aes(Timestamp, residuals_up_down), 
             color = "#000000",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp > max(calib_time) &
                        Timestamp %in% TN_time),
             aes(Timestamp, residuals_up_down),
             color = "#999999",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% TP_time),
             aes(Timestamp, residuals_up_down),
             color = "#009E73",
             shape = 17,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FP_time),
             aes(Timestamp, residuals_up_down),
             color = "#E69F00",
             shape = 18,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FN_time),
             aes(Timestamp, residuals_up_down),
             color = "#FF3333",
             shape = 15,
             size = 2) +
  ylab("residuals")


###########################
  ## Fitting-GAM-down-AR
###########################

#adding temperature and conductance at lag1 from the downstream sensor 
#to the GAM-up

gam_down_AR <- gam(turbidity_downstream_log ~
                  s(turbidity_downstream_log_lag1) +
                  s(conductance_downstream) +
                  s(dissolved_oxygen_downstream) +
                  s(level_downstream, k=6) +
                  s(temperature_downstream),
                data = data_train)




## SPOT-for-GAM-down-AR

# computing residuals using the model 3 from all data

data_residuals <- data_residuals %>%
  mutate(predict_mod_down_AR = as.numeric(predict.gam(gam_down_AR, newdata = data_all)),
         residuals_down_AR = turbidity_downstream_log - predict_mod_down_AR)


confusion_mod_down_AR <- compute_confusion(data = data_residuals, 
                                           residuals = "residuals_down_AR",
                                           name = "model_down_AR", n = n,
                                           q = q, t_prob = t_prob,
                                           cmax = FALSE)

confusion_mod_down_AR <- confusion_mod_down_AR %>%
  rename(AnomalyFlag_SPOT_down_AR = AnomalyFlag_SPOT,
         confusion_mod_down_AR = confusion_matrix)
data_residuals <- data_residuals %>%
  left_join(confusion_mod_down_AR, by = "Timestamp")



## classification-GAM-down-AR

classif_plot_turb_down_AR <- data_residuals %>%
  select(Timestamp, turbidity_downstream, confusion_mod_down_AR) %>%
  mutate(confusion_mod_down_AR = factor(confusion_mod_down_AR,
                                        levels = c("TP", "TN", "FP", "FN",
                                                   "NA"))) %>%
  drop_na() %>%
  ggplot() +
  geom_point(aes(Timestamp, turbidity_downstream,
                 color = confusion_mod_down_AR,
                 shape = confusion_mod_down_AR,
                 size = confusion_mod_down_AR,
                 alpha = confusion_mod_down_AR)) +
  scale_color_manual(values = c("#009E73", "#999999", "#E69F00",
                                "#FF3333")) +
  scale_shape_manual(values = c(17,16,18,15)) +
  scale_size_manual(values = c(2,0.3,2,2)) +
  scale_alpha_manual(values = c(1,1,1,1)) +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  ylab("Turbidity (FNU)")

calib_time <- data_residuals %>% 
  slice(1:n) %>% 
  pull(Timestamp)

TP_time <- data_residuals %>% 
  filter(confusion_mod_down_AR == "TP") %>% 
  pull(Timestamp)

TN_time <- data_residuals %>% 
  filter(confusion_mod_down_AR == "TN") %>% 
  pull(Timestamp)

FP_time <- data_residuals %>% 
  filter(confusion_mod_down_AR == "FP") %>% 
  pull(Timestamp)

FN_time <- data_residuals %>% 
  filter(confusion_mod_down_AR == "FN") %>% 
  pull(Timestamp)


classif_plot_resid_down_AR <- data_residuals %>%
  select(Timestamp, residuals_down_AR) %>%
  drop_na() %>%
  ggplot() +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% calib_time),
             aes(Timestamp, residuals_down_AR), 
             color = "#000000",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp > max(calib_time) &
                        Timestamp %in% TN_time),
             aes(Timestamp, residuals_down_AR),
             color = "#999999",
             shape = 16,
             size = 0.3,
             alpha = 0.3) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% TP_time),
             aes(Timestamp, residuals_down_AR),
             color = "#009E73",
             shape = 17,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FP_time),
             aes(Timestamp, residuals_down_AR),
             color = "#E69F00",
             shape = 18,
             size = 2) +
  geom_point(data = data_residuals %>% 
               filter(Timestamp %in% FN_time),
             aes(Timestamp, residuals_down_AR),
             color = "#FF3333",
             shape = 15,
             size = 2) +
  ylab("residuals")




## ---- classification-all-models

plot_turb_classif <- data_residuals %>% 
  select(Timestamp, turbidity_downstream, turbidityAnomalyFlag) %>% 
  ggplot(aes(Timestamp, turbidity_downstream)) +
  geom_point(size = 0.3,
             color = "#999999",
             shape = 16) +
  geom_point(data = data_residuals %>% 
               filter(turbidityAnomalyFlag == "outlier"),
             aes(Timestamp, turbidity_downstream),
             color = "red",
             size = 2,
             shape = 17) +
  ylab("Turbidity")

grid.arrange(arrangeGrob(plot_turb_classif,
                         classif_plot_resid_up +
                           ylab("GAM-up"),
                         classif_plot_resid_up_AR +
                           ylab("GAM-up-AR"),
                         classif_plot_resid_down +
                           ylab("GAM-down"),
                         classif_plot_resid_down_AR +
                           ylab("GAM-down-AR"),
                         classif_plot_resid_up_down +
                           ylab("GAM-up-down"),
                         ncol = 1),
             ncol = 1, mylegend, heights=c(10, 1))



## ---- model-comparison

confussion_matrix <- data_residuals %>%
  drop_na(turbidity_downstream_log) %>%
  dplyr::select(Timestamp, confusion_mod_up_AR,
                confusion_mod_up, confusion_mod_down,
                confusion_mod_up_down, 
                confusion_mod_down_AR) %>%
  pivot_longer(-Timestamp, names_to = "method", values_to = "value") %>%
  group_by(method) %>%
  count(value) %>%
  spread(key = value, value = n) %>%
  mutate(method = recode(method,
                         confusion_mod_up_AR = "GAM-up-AR",
                         confusion_mod_up = "GAM-up",
                         confusion_mod_down = "GAM-down",
                         confusion_mod_up_down = "GAM-up-down",
                         confusion_mod_down_AR = "GAM-down-AR"))


tbl_comparison <- confussion_matrix %>%
  mutate(FN = ifelse(is.na(FN), 0, FN),
         FP = ifelse(is.na(FP), 0, FP),
         accuracy = round((TP + TN)/(TP + TN + FP + FN), digits = 4),
         ER = round((FP + FN)/(TP + TN + FP + FN), digits = 4),
         NPV = round((TN)/(TN + FN), digits = 4),
         PPV = round((TP)/(TP + FP), digits = 4),
         Sp = TN/(TN+FP),
         Sn = TP/(TP+FN),
         Nn = (FP+TN)/(TP + TN + FP + FN),
         Np = (TP + FN)/(TP + TN + FP + FN),
         P = Sp*Nn + Sn*Np,
         RI = abs(Sp-Sn)/(Sp+Sn),
         OP = round(P - RI, digits = 4)) %>%
  select(method, TP, TN, FP, FN, accuracy, ER, OP, NPV, PPV) %>%
  arrange(desc(OP))

params <- data.frame(
  method = c("GAM-up-AR", "GAM-up", "GAM-down", "GAM-up-down", 
             "GAM-down-AR"),
  Predictors = c("turbidity(up), conductance(up), level(up),
                 temperature(up), lagged responses",
                 "turbidity(up), conductance(up), level(up),
                 temperature(up), dissolved oxygen(up)",
                 "conductance(down), level(down), temperature(down),
                 dissolved oxygen(down)",
                 "turbidity(up), conductance(up), level(up),
                 temperature(up),dissolved oxygen(up),conductance(down),
                 temperature(down)",
                 "conductance(down), dissolved oxygen(down),
                 level(down), temperature(down), lagged response"
  ))

options(scipen = 999)

tbl_comparison <- tbl_comparison %>%
  left_join(params, by = "method") %>%
  select(method, Predictors, everything())



########################################
    ##-- oddwater comparison --##
########################################

## ---- oddwater-comparison

PRIN_5min_turb <- PRIN_5min_flagged %>%
  filter(site == "down") %>%
  select(roundedTimestamp, turbidity, turbidityAnomalyFlag  ) %>%
  rename(
    Timestamp = roundedTimestamp,
    flag = turbidityAnomalyFlag
  ) %>%
  mutate(Timestamp = ymd_hms(Timestamp),
         flag = as_factor(flag)) %>%
  drop_na()

turb_oct <- PRIN_5min_turb %>%
  filter(Timestamp >= ymd("2019-10-01") & 
           Timestamp < ymd("2020-01-01"))

# A function to do different transformations
transform_data1 <- function(data,  time_bound = 90, regular = FALSE, time_col = "Timestamp")
{
  if(!(lubridate::is.Date(data[[time_col]]) | lubridate::is.POSIXct(data[[time_col]]) ))
  { data[time_col] <- lubridate::dmy_hm((data[[time_col]]) )}
  
  n <- nrow(data)
  data_var <- as.matrix(data[ , !(names(data) %in% time_col)])
  
  # apply log transformation
  log_series <- log(data_var)
  colnames(log_series) <- paste("log_", colnames(log_series), sep = "")
  data <- cbind(data, log_series)
  
  # take the first difference of the original series
  diff_series <- rbind( rep(NA, ncol(data_var)), diff(data_var))
  colnames(diff_series) <- paste("diff_", colnames(diff_series), sep = "")
  data <- cbind(data, diff_series)
  
  # take the first difference of the log series
  diff_log_series <- rbind( rep(NA, ncol(data_var)), diff(log_series))
  colnames(diff_log_series) <- paste("difflog_", colnames(diff_log_series), sep = "")
  data <- cbind(data, diff_log_series)
  
  
  # rate of change
  rc_series <- (data_var[2:n,] - data_var[1:(n-1),]) / data_var[1:(n-1),] %>% as.matrix(col=1)
  rc_series <-  rbind( NA,rc_series)
  colnames(rc_series) <- paste("rc_", colnames(data_var), sep = "")
  data <- cbind(data, rc_series)
  
  # Ratio
  ratio_series <-  (data_var[2:n,] / data_var[1:(n-1),]) %>%
    as.matrix(col=1)
  ratio_series <-  rbind( NA,ratio_series)
  colnames(ratio_series) <- paste("ratio_", colnames(data_var), sep = "")
  data <- cbind(data, ratio_series)
  
  # Relative difference (log)
  relative_series <-  (log_series[2:(n-1),] - (1/2)*(log_series[1:(n-2),] +log_series[3:n,]  )) %>% as.matrix(col=1)
  
  relative_series <-  rbind( NA,relative_series, NA)
  colnames(relative_series) <- paste("rdifflog_", colnames(data_var), sep = "")
  data <- cbind(data, relative_series)
  
  # Relative difference (original)
  relative_series_o <- (data_var[2:(n-1),] - (1/2)*(data_var[1:(n-2),] +data_var[3:n,]  ))%>% as.matrix(col=1)
  relative_series_o <-  rbind( NA,relative_series_o, NA)
  colnames(relative_series_o) <- paste("rdiff_", colnames(data_var), sep = "")
  data <- cbind(data, relative_series_o)
  
  data <- tsibble::as_tsibble(data, index = time_col , regular = regular)
  
  return(data)
}


data <- turb_oct %>%
  select(-flag) %>%
  tidyr::drop_na() 
trans_data <-transform_data1(data)
full_data <- left_join(trans_data, turb_oct %>% select(-turbidity), by = "Timestamp")

code <- function(x, y) {
  if (x == "1" & y == "outlier") {
    return("TP")
  }
  if (x == "1" & y == "typical") {
    return("FN")
  }
  if (x == "0" & y == "outlier") {
    return("FP")
  }
  if (x == "0" & y == "typical") {
    return("TN")
  }
}

find_performance<- function(data, p, alpha ,k)
{
  out_rc <- stray::find_HDoutliers(data[["rc_turbidity"]], p= p,
                                   alpha = alpha, k=k)
  
  turb_data <- data %>%
    select(Timestamp, turbidity, flag) %>%
    mutate(type = out_rc$type)
  
  code_trub <- mapply(
    code, turb_data[["flag"]],
    turb_data[["type"]]
  )  
  
  turb_data <- turb_data %>%
    mutate(code  = code_trub)
  
  
  sum_tub <- list(p=p, alpha = alpha, k=k, 
                  FN = sum(turb_data$code=="FN") ,  
                  FP = sum(turb_data$code=="FP") ,  
                  TN = sum(turb_data$code=="TN"),  
                  TP = sum(turb_data$code=="TP"))
  
  return(sum_tub %>% rbind() %>% as_tibble())
  
  
}


sumtab <- find_performance(data =full_data, p=0.5, alpha=0.01 ,k=1)

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.5, alpha=0.01 ,k=5))

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.75, alpha=0.01 ,k=1))

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.75, alpha=0.01 ,k=5))


stray_spec <- data.frame(
  method = c("stray(p=0.5, k=1)", "stray(p=0.5, k=5)", "stray(p=0.75, k=1)",
             "stray(p=0.75, k=5)"),
  Predictors = c("p=0.5, k=1", "p=0.5, k=5", "p=0.75, k=1", "p=0.75, k=5"))

sumtab <- sumtab %>% 
  select(TP, TN, FP, FN) %>% 
  bind_cols(stray_spec) %>% 
  select(method, Predictors, everything()) %>% 
  mutate_at(.vars = c("TP", "TN", "FP", "FN"), 
            .funs = as.numeric) %>% 
  mutate(FN = ifelse(is.na(FN), 0, FN),
         FP = ifelse(is.na(FP), 0, FP),
         accuracy = round((TP + TN)/(TP + TN + FP + FN), digits = 4),
         ER = round((FP + FN)/(TP + TN + FP + FN), digits = 4),
         NPV = round((TN)/(TN + FN), digits = 4),
         PPV = round((TP)/(TP + FP), digits = 4),
         Sp = TN/(TN+FP),
         Sn = TP/(TP+FN),
         Nn = (FP+TN)/(TP + TN + FP + FN),
         Np = (TP + FN)/(TP + TN + FP + FN),
         P = Sp*Nn + Sn*Np,
         RI = abs(Sp-Sn)/(Sp+Sn),
         OP = round(P - RI, digits = 4)) %>%
  select(method, Predictors, TP, TN, FP, FN, accuracy, ER, OP, NPV, PPV)

# print(sumtab)


## ---- all-comparison

# tbl_comparison <- tbl_comparison %>% 
#   bind_rows(sumtab) %>% 
#   arrange(desc(OP)) %>% 
#   select(-Predictors)

# print(tbl_comparison)

# tbl_comparison <- tbl_comparison %>%
#   bind_rows(sumtab) %>%
#   arrange(desc(OP)) %>% 
#   dplyr::select(method, TP, TN, FP, FN, 
#                 `OP`)
# 
# tbl_comparison
# save(tbl_comparison, file = here::here("ARCLP-DES-Catchup", "Meeting-4",
#                                        "plots", "tbl_comparison.rda"))
# 


tbl_comparison %>%
  bind_rows(sumtab) %>%
  arrange(desc(OP)) %>%
  dplyr::select(method, TP, TN, FP, FN,
                "OP") %>%
  kableExtra::kbl(format = "latex", booktabs = T,
                    linesep = "") %>%
  kableExtra::kable_styling(font_size = 7,
                            latex_options="scale_down") 
# %>% 
#   kableExtra::column_spec(1, bold = TRUE) %>% 
#   kableExtra::row_spec(0, bold = TRUE) 
  



#%>%
#   kableExtra::column_spec(1, bold = T) %>%
#   kableExtra::row_spec(0, bold = T)

# kableExtra::kbl(tbl_comparison, booktabs = TRUE,
#                 format = "latex",
#                 linesep = "\\addlinespace") %>%
#   kableExtra::kable_styling(latex_options="scale_down")
#   # kableExtra::column_spec(2, width = "4cm")
#   # kableExtra::column_spec(c(1,2,8), width = c("1.5cm","4cm","1cm"))




