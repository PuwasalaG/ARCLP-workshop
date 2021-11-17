## ----loadLibraries
library(tidyverse)
library(lubridate)
library(mgcv)

# ----loadFunctions
source(here::here("ARCLP-DES-Catchup", "Meeting-3", "scripts", "SPOT_algorithm.R")
)

# ----loadData
load(here::here("ARCLP-DES-Catchup", "Meeting-3", "data", "data_upstream_lagged.rda")
)
load(here::here("ARCLP-DES-Catchup", "Meeting-3", "data", "data_downstream.rda")
)
load(here::here("ARCLP-DES-Catchup", "Meeting-3", "data", "PRIN_5min_cleaned.rda")
)

## ----PrepareData
data_downstream <- data_downstream %>% 
  rename_at(.vars = vars(-c(Timestamp, turbidityAnomalyFlag)),
            .funs = ~paste(., "_downstream", sep = ""))

#Since the turbidity is too skewed we take the log of turbidity
#and computing the lags of downstream variables
data <- data_upstream_lagged %>% 
  select(-turbidity_downstream, -turbidityAnomalyFlag) %>% 
  left_join(data_downstream, by = "Timestamp") %>% 
  mutate(turbidity_downstream_log = log(turbidity_downstream),
         turbidity_upstream_dt_log = log(turbidity_upstream_dt),
         level_downstream_lag1 = lag(level_downstream, n = 1L),
         temperature_downstream_lag1 = lag(temperature_downstream, 
                                           n = 1L),
         conductance_downstream_lag1 = lag(conductance_downstream, 
                                           n = 1L),
         turbidityAnomalyFlag = factor(turbidityAnomalyFlag,
                                       levels = c(0,1),
                                       labels = c("typical",
                                                  "outlier")))

## ----TimePlots
PRIN_5min_cleaned %>% 
  select(roundedTimestamp, specificConductance, turbidity, 
         surfacewaterElevMean, surfWaterTempMean, site) %>% 
  rename("Timestamp" = roundedTimestamp,
         "Conductance" = specificConductance,
         "Turbidity" = turbidity,
         "Level" = surfacewaterElevMean,
         "Temperature" = surfWaterTempMean) %>%
  filter(Timestamp >= ymd("2019-10-01") & Timestamp < ymd("2020-01-01")) %>% 
  pivot_longer(-c(Timestamp, site)) %>% 
  ggplot(aes(Timestamp, value, color = site)) +
  geom_line() +
  facet_wrap(~name, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("#999999", "#E69F00"))


data %>% 
  select(Timestamp, turbidity_downstream, turbidityAnomalyFlag) %>% 
  ggplot(aes(Timestamp, turbidity_downstream)) +
  geom_point() +
  geom_point(data = data %>% 
               filter(turbidityAnomalyFlag == "outlier"),
             aes(Timestamp, turbidity_downstream, color = "red")) +
  theme(legend.position = "none")


## ----trainingData

data_train <- data %>% 
  filter(Timestamp >= ymd("2019-10-01") &
           Timestamp < ymd("2019-12-01"))

# We will remove the anomalies in the turbidity in training data
data_train <- data_train %>% 
  mutate(turbidity_downstream = if_else(turbidityAnomalyFlag=="outlier", as.numeric(NA), 
                                        turbidity_downstream))


## ----Model1

Model1 <- gam(turbidity_downstream_log ~ 
                s(turbidity_upstream_dt_log) +
                s(level_upstream_dt) +
                s(conductance_upstream_dt) +
                s(temperature_upstream_dt),
              data = data_train)

#computing residuals from the model1 for all data
data <- data %>% 
  mutate(pred_turbidity_downstream_log_model1 =
           as.numeric(predict.gam(Model1,
                                  newdata = data)),
         residuals_model1 = turbidity_downstream_log -
           pred_turbidity_downstream_log_model1)


## ----Model2

Model2 <- gam(turbidity_downstream_log ~ 
                s(level_downstream_lag1) +
                s(conductance_downstream_lag1) +
                s(temperature_downstream_lag1),
              data = data_train)

#computing residuals from the model2
data <- data %>% 
  mutate(pred_turbidity_downstream_log_model2 =
           as.numeric(predict.gam(Model2,
                                  newdata = data)),
         residuals_model2 = turbidity_downstream_log -
           pred_turbidity_downstream_log_model2)


## ----Model3

Model3 <- update(Model1, .~.+s(level_downstream_lag1)+
                   s(conductance_downstream_lag1) +
                   s(temperature_downstream_lag1))


#computing residuals from the model3 for all data
data <- data %>% 
  mutate(pred_turbidity_downstream_log_model3 =
           as.numeric(predict.gam(Model3,
                                  newdata = data)),
         residuals_model3 = turbidity_downstream_log -
           pred_turbidity_downstream_log_model3)


data %>% 
  select(Timestamp, turbidity_downstream, residuals_model1,
         residuals_model2, residuals_model3) %>% 
  pivot_longer(-Timestamp) %>%
  mutate(name = factor(name, levels = c("turbidity_downstream", 
                                        "residuals_model1",
                                        "residuals_model2", 
                                        "residuals_model3"))) %>% 
  ggplot(aes(Timestamp, value)) +
  geom_point(size = 0.5) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scales = "free_y")


## ----outlierThresholdComputation

q <- 0.001
t_prob <- 0.98
n = 24*24*60/5 # 24 days for calibration

# from model 1
residuals_model1 <- data %>% 
  select(Timestamp, residuals_model1) %>% 
  rename("residuals" = residuals_model1) %>% 
  mutate(residuals = as.numeric(residuals))

residuals_model1 <- find_outliers_SPOT(residuals = residuals_model1,
                                       n = n,
                                       initial_threshold_prob = t_prob,
                                       level = q)

residuals_model1 <- residuals_model1 %>% 
  select(-residuals) %>% 
  rename("residuals_model1_anomaly" = residuals_anomaly_SPOT)


# from model 2
residuals_model2 <- data %>% 
  select(Timestamp, residuals_model2) %>% 
  rename("residuals" = residuals_model2) %>% 
  mutate(residuals = as.numeric(residuals))

residuals_model2 <- find_outliers_SPOT(residuals = residuals_model2,
                                       n = n,
                                       initial_threshold_prob = t_prob,
                                       level = q)

residuals_model2 <- residuals_model2 %>% 
  select(-residuals) %>% 
  rename("residuals_model2_anomaly" = residuals_anomaly_SPOT)


# from model 3
residuals_model3 <- data %>% 
  select(Timestamp, residuals_model3) %>% 
  rename("residuals" = residuals_model3) %>% 
  mutate(residuals = as.numeric(residuals))

residuals_model3 <- find_outliers_SPOT(residuals = residuals_model3,
                                       n = n,
                                       initial_threshold_prob =
                                         t_prob,
                                       level = q)

residuals_model3 <- residuals_model3 %>% 
  select(-residuals) %>% 
  rename("residuals_model3_anomaly" = residuals_anomaly_SPOT)


data <- list(data, residuals_model1, residuals_model2,
             residuals_model3) %>%
  purrr::reduce(left_join, by = "Timestamp")


## ----confussionMatrix


data <- data %>% 
  mutate(confusion_matrix_model1 = case_when(turbidityAnomalyFlag
                                             == "typical" 
                                             & residuals_model1_anomaly ==
                                               "outlier" ~ "FP",
                                             turbidityAnomalyFlag
                                             == "outlier" 
                                             & residuals_model1_anomaly ==
                                               "typical" ~ "FN",
                                             turbidityAnomalyFlag
                                             == "outlier" 
                                             & residuals_model1_anomaly ==
                                               "outlier" ~ "TP",
                                             turbidityAnomalyFlag
                                             == "typical" 
                                             & residuals_model1_anomaly ==
                                               "typical" ~ "TN",
                                             TRUE ~
                                               as.character(NA)),
         confusion_matrix_model1 = factor(confusion_matrix_model1, 
                                          levels = c("TP", "TN", "FP", "FN",
                                                     "NA")),
         confusion_matrix_model2 = case_when(turbidityAnomalyFlag
                                             == "typical" 
                                             & residuals_model2_anomaly ==
                                               "outlier" ~ "FP",
                                             turbidityAnomalyFlag
                                             == "outlier" 
                                             & residuals_model2_anomaly ==
                                               "typical" ~ "FN",
                                             turbidityAnomalyFlag
                                             == "outlier" 
                                             & residuals_model2_anomaly ==
                                               "outlier" ~ "TP",
                                             turbidityAnomalyFlag
                                             == "typical" 
                                             & residuals_model2_anomaly ==
                                               "typical" ~ "TN",
                                             TRUE ~
                                               as.character(NA)),
         confusion_matrix_model2 = factor(confusion_matrix_model2, 
                                          levels = c("TP", "TN", "FP", "FN",
                                                     "NA")),
         confusion_matrix_model3 = case_when(turbidityAnomalyFlag
                                             == "typical" 
                                             & residuals_model3_anomaly ==
                                               "outlier" ~ "FP",
                                             turbidityAnomalyFlag
                                             == "outlier" 
                                             & residuals_model3_anomaly ==
                                               "typical" ~ "FN",
                                             turbidityAnomalyFlag
                                             == "outlier" 
                                             & residuals_model3_anomaly ==
                                               "outlier" ~ "TP",
                                             turbidityAnomalyFlag
                                             == "typical" 
                                             & residuals_model3_anomaly ==
                                               "typical" ~ "TN",
                                             TRUE ~
                                               as.character(NA)),
         confusion_matrix_model3 = factor(confusion_matrix_model3, 
                                          levels = c("TP", "TN", "FP", "FN",
                                                     "NA"))
  )


## ----Comparison

confussion_matrix <- data %>% 
  drop_na(turbidity_downstream_log) %>% 
  dplyr::select(Timestamp, confusion_matrix_model1, 
                confusion_matrix_model2,
                confusion_matrix_model3) %>% 
  pivot_longer(-Timestamp, names_to = "method", values_to = "value") %>% 
  group_by(method) %>% 
  count(value) %>% 
  spread(key = value, value = n) %>% 
  mutate(method = recode(method, 
                         confusion_matrix_model1 = "model1",
                         confusion_matrix_model2 = "model2",
                         confusion_matrix_model3 = "model3"))

confussion_matrix <- confussion_matrix %>% 
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
  select(TP, TN, FP, FN, OP, accuracy, ER) %>% 
  arrange(desc(OP))  

data %>% 
  select(Timestamp, turbidity_downstream, 
         confusion_matrix_model3) %>% 
  drop_na() %>% 
  ggplot() +
  geom_point(aes(Timestamp, turbidity_downstream, 
                 color = confusion_matrix_model3,
                 shape = confusion_matrix_model3)) +
  scale_color_manual(values = c("#009E73", "#999999", "#E69F00",
                                "#FF3333","black")) +
  scale_shape_manual(values = c(17,16,18,15,16)) +
  scale_size_manual(values = c(2,0.3,2,2,0.3)) +
  theme(strip.text = element_text(size = 8)) +
  ggtitle("Classifying the outliers detected from model−3")



save(confussion_matrix, file = here::here("ARCLP-DES-Catchup", "Meeting-3", "plots",
                                     "confussion_matrix.rda"))


