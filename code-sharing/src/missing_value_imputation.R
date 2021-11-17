############################################
    ##-- Missing value imputation --##
############################################

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

## ---- loadPRINData
# Pringle Creek in Texas
load(here::here("code-sharing", "data", "PRIN_5min_cleaned.rda"))
PRIN_5min_cleaned <- PRIN_5min_cleaned %>%
  rename(
    Timestamp = roundedTimestamp,
    level = surfacewaterElevMean,
    conductance = specificConductance,
    dissolved_oxygen = dissolvedOxygen,
    temperature = surfWaterTempMean
  ) %>%
  mutate(Timestamp = ymd_hms(Timestamp))


## ---- dataPreProccessing

# Choosing data from "2019-10-01" to "2019-12-31"
PRIN_5min_new <- PRIN_5min_cleaned %>%
  filter(Timestamp >= ymd("2019-10-01") 
         & Timestamp < ymd("2020-01-01")) %>% 
  select(str_subset(names(PRIN_5min_cleaned), "Anomaly", negate = T))

# selecting turbidity from the downstream and level, temperature, conductance
# dissolved_oxygen and pH from the upstream sensor
# from upstream sensor
data_upstream <- PRIN_5min_new %>% 
  filter(site == "up") %>% 
  select(Timestamp, level, temperature, conductance, dissolved_oxygen,
         pH, turbidity) %>% 
  rename("level_upstream" = level,
         "temperature_upstream" = temperature,
         "conductance_upstream" = conductance,
         "dissolved_oxygen_upstream" = dissolved_oxygen,
         "pH_upstream" = pH,
         "turbidity_upstream" = turbidity) %>% 
  mutate(turbidity_upstream = if_else(turbidity_upstream == 0, 
                                      as.numeric(NA), 
                                      turbidity_upstream))

data_downstream <- PRIN_5min_new %>% 
  filter(site == "down") %>% 
  select(Timestamp, level, temperature, conductance, dissolved_oxygen,
         pH, turbidity) %>% 
  rename("level_downstream" = level,
         "temperature_downstream" = temperature,
         "conductance_downstream" = conductance,
         "dissolved_oxygen_downstream" = dissolved_oxygen,
         "pH_downstream" = pH,
         "turbidity_downstream" = turbidity) %>% 
  mutate(turbidity_downstream = if_else(turbidity_downstream == 0, 
                                        as.numeric(NA), 
                                        turbidity_downstream),
         conductance_downstream = if_else(Timestamp == ymd_hms("2019-12-09 15:40:00"),
                                          as.numeric(NA), conductance_downstream))


## ---- missing_value_imputation 

# We impute the missing values in the upstream variables
## --Imputation - Upstream sensor--

#level - using linear interpolation
#temperature, conductance, pH and DO - using Kalman filter with ARIMA state-space models
#turbidity - conditional normalisation followed by kalman filter with ARIMA

impute_predicts <- function(x){
  fit <- auto.arima(x, stepwise = F)
  model <- fit$model
  x_impute <- na_kalman(x, model = model)
  return(x_impute)
}

data_upstream_imputed <- data_upstream %>% 
  mutate_at(.vars = vars(temperature_upstream, conductance_upstream,
                         dissolved_oxygen_upstream, pH_upstream),
            .funs = ~impute_predicts(.x)) %>% 
  mutate(level_upstream = na_interpolation(level_upstream))


#Imputing missing values in turbidity upstream using conditional normalisation methods
# refer "conditional normalisation in time series analysis" paper

# first computing conditional moments
# We use only level and temperature as predictors from the upstream sensor
cond_moments_turb_up <- conditional_moments(data = data_upstream_imputed,
                                            x = turbidity_upstream,
                                            z_numeric =
                                              c(level_upstream,
                                                temperature_upstream),
                                            knots_mean = c(5,5),
                                            knots_variance =
                                              c(4,4))

autoplot(cond_moments_turb_up, type = "mean")
autoplot(cond_moments_turb_up, type = "variance")
summary(cond_moments_turb_up, type = "mean")
summary(cond_moments_turb_up, type = "variance")

# conditional normalisation
data_norm_turb_up <- cond_moments_turb_up$data_conditional_moments

data_norm_turb_up <- data_norm_turb_up %>%
  mutate(normalised_turbidity_upstream =
           as.numeric((turbidity_upstream - E_turbidity_upstream)/sqrt(Var_turbidity_upstream)))


# imputing missing values in normalised turbidity with a Kalman smoother
# here we choose the arima model using auto.arima
x <- data_norm_turb_up$normalised_turbidity_upstream

fit_arima <- auto.arima(x, stepwise = F)
model <- fit_arima$model
x_impute <- na_kalman(x, model = model)

# ggplot_na_imputations(x, x_impute, size_points = 0.5, 
#                       title = "Kalman smoothing")

data_norm_turb_up <- data_norm_turb_up %>% 
  mutate(normalised_turbidity_upstream_impute = x_impute,
         turbidity_upstream_imputed =
           if_else(is.na(turbidity_upstream), normalised_turbidity_upstream_impute*sqrt(Var_turbidity_upstream)+E_turbidity_upstream, turbidity_upstream))

# data_norm_turb_up %>% 
#   select(Timestamp, turbidity_upstream, turbidity_upstream_imputed) %>% 
#   pivot_longer(-Timestamp) %>% 
#   ggplot(aes(Timestamp, value)) +
#   geom_line() +
#   facet_wrap(~name, ncol = 1)

data_upstream_imputed <- data_upstream_imputed %>% 
  mutate(turbidity_upstream = data_norm_turb_up$turbidity_upstream_imputed,
         turbidity_upstream = if_else(turbidity_upstream<=0, 
                                      as.numeric(NA), 
                                      turbidity_upstream))


save(data_upstream_imputed, file = here::here("code-sharing",
                                              "data",
                                              "data_upstream_imputed.rda"))

## --Imputation - Downstream sensor--

#level - using linear interpolation
#temperature, conductance, pH and DO - using Kalman filter with ARIMA state-space models
#turbidity - conditional normalisation followed by kalman filter with ARIMA

data_downstream_imputed <- data_downstream %>%
  mutate_at(.vars = vars(temperature_downstream, conductance_downstream,
                         dissolved_oxygen_downstream, pH_downstream),
            .funs = ~impute_predicts(.x)) %>%
  mutate(level_downstream = na_interpolation(level_downstream))


#Imputing missing values in turbidity downstream using conditional normalisation methods
# refer "conditional normalisation in time series analysis" paper

# first computing conditional moments
# We use only level and temperature as predictors from the upstream sensor
cond_moments_turb_down <- conditional_moments(data = data_downstream_imputed,
                                              x = turbidity_downstream,
                                              z_numeric =
                                                c(level_downstream,
                                                  temperature_downstream),
                                              knots_mean = c(5,5),
                                              knots_variance =
                                                c(4,4))

autoplot(cond_moments_turb_down, type = "mean")
autoplot(cond_moments_turb_down, type = "variance")
summary(cond_moments_turb_down, type = "mean")
summary(cond_moments_turb_down, type = "variance")

# conditional normalisation
data_norm_turb_down <- cond_moments_turb_down$data_conditional_moments

data_norm_turb_down <- data_norm_turb_down %>%
  mutate(normalised_turbidity_downstream =
           as.numeric((turbidity_downstream - E_turbidity_downstream)/sqrt(Var_turbidity_downstream)))


# imputing missing values in normalised turbidity with a Kalman smoother
# here we choose the arima model using auto.arima
y <- data_norm_turb_down$normalised_turbidity_downstream

fit_arima <- auto.arima(y, stepwise = F)
model <- fit_arima$model
y_impute <- na_kalman(y, model = model)

# ggplot_na_imputations(y, y_impute, size_points = 0.5,
#                       title = "Kalman smoothing")

data_norm_turb_down <- data_norm_turb_down %>%
  mutate(normalised_turbidity_downstream_impute = y_impute,
         turbidity_downstream_imputed =
           if_else(is.na(turbidity_downstream), normalised_turbidity_downstream_impute*sqrt(Var_turbidity_downstream)+E_turbidity_downstream, turbidity_downstream))

# data_norm_turb_down %>%
#   select(Timestamp, turbidity_downstream, turbidity_downstream_imputed) %>%
#   pivot_longer(-Timestamp) %>%
#   ggplot(aes(Timestamp, value)) +
#   geom_line() +
#   facet_wrap(~name, ncol = 1)

data_downstream_imputed <- data_downstream_imputed %>%
  mutate(turbidity_downstream = data_norm_turb_down$turbidity_downstream_imputed,
         turbidity_downstream = if_else(turbidity_downstream<=0,
                                        as.numeric(NA),
                                        turbidity_downstream))

