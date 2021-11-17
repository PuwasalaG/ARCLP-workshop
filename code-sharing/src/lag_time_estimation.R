### Supplemental Materials - 1 (Data and scripts)
### Lag time estimation and computing lagged upstream data

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
load(here::here("weird_river_data_paper", "data", "PRIN_5min_cleaned.rda"))
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




## ---- timeplots-upstream

is_na <- function(x){
  x_is.na = if_else(is.na(x)|is.nan(x), 1, 0)
  x_is.na = factor(x_is.na, levels = c(0,1))
  return(x_is.na)
}

data_upstream_is.na <- data_upstream %>%
  mutate(across(.cols = c(level_upstream, temperature_upstream, 
                          conductance_upstream, 
                          dissolved_oxygen_upstream, pH_upstream, 
                          turbidity_upstream),
                .fns = ~is_na(.x)))

p_level_up <- data_upstream_imputed %>% 
  select(Timestamp, level_upstream) %>% 
  left_join(y = data_upstream_is.na %>% 
              select(Timestamp, level_upstream), 
            by = "Timestamp") %>% 
  ggplot(aes(Timestamp, level_upstream.x, color = level_upstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_upstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, level_upstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"), 
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("Level (m)") +
  ggtitle("Level - upstream")

p_temp_up <- data_upstream_imputed %>% 
  select(Timestamp, temperature_upstream) %>% 
  left_join(y = data_upstream_is.na %>% 
              select(Timestamp, temperature_upstream), 
            by = "Timestamp") %>% 
  ggplot(aes(Timestamp, temperature_upstream.x, 
             color = temperature_upstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_upstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, temperature_upstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"), 
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab(expression('Temperature ('~degree*C*')')) +
  ggtitle("Temperature - upstream")

p_cond_up <- data_upstream_imputed %>% 
  select(Timestamp, conductance_upstream) %>% 
  left_join(y = data_upstream_is.na %>% 
              select(Timestamp, conductance_upstream), 
            by = "Timestamp") %>% 
  ggplot(aes(Timestamp, conductance_upstream.x, 
             color = conductance_upstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_upstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, conductance_upstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"), 
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab(expression(Conductance~(mu*S*'/'*cm))) +
  ggtitle("Conductance - upstream")

p_turb_up <- data_upstream_imputed %>% 
  select(Timestamp, turbidity_upstream) %>% 
  left_join(y = data_upstream_is.na %>% 
              select(Timestamp, turbidity_upstream), 
            by = "Timestamp") %>% 
  ggplot(aes(Timestamp, turbidity_upstream.x, 
             color = turbidity_upstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_upstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, turbidity_upstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"), 
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("Turbidity (FNU)") +
  ggtitle("Turbidity - upstream")

p_DO_up <- data_upstream_imputed %>% 
  select(Timestamp, dissolved_oxygen_upstream) %>% 
  left_join(y = data_upstream_is.na %>% 
              select(Timestamp, dissolved_oxygen_upstream), 
            by = "Timestamp") %>% 
  ggplot(aes(Timestamp, dissolved_oxygen_upstream.x, 
             color = dissolved_oxygen_upstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_upstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, dissolved_oxygen_upstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"), 
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("DO (mg/l)") +
  ggtitle("Dissolved Oxygen - upstream")

p_pH_up <- data_upstream_imputed %>% 
  select(Timestamp, pH_upstream) %>% 
  left_join(y = data_upstream_is.na %>% 
              select(Timestamp, pH_upstream), 
            by = "Timestamp") %>% 
  ggplot(aes(Timestamp, pH_upstream.x, 
             color = pH_upstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_upstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, pH_upstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"), 
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("pH") +
  ggtitle("pH - upstream")

# p_turb_down <- data_downstream %>% 
#   select(Timestamp, turbidity_downstream) %>% 
#   ggplot(aes(Timestamp, turbidity_downstream)) +
#   geom_point(size = 0.5) +
#   scale_colour_manual(values = c("#000000", "#D55E00"), 
#                       breaks = c(0,1)) +
#   theme(legend.position = "none") +
#   ylab("Turbidity (FNU)") +
#   ggtitle("Turbidity - downstream")

ggpubr::ggarrange(p_turb_up, p_cond_up, p_DO_up, p_pH_up, p_level_up,
                  p_temp_up, ncol = 1)





## ---- timeplots-downstream

data_downstream_is.na <- data_downstream %>%
  mutate(across(.cols = c(level_downstream, temperature_downstream,
                          conductance_downstream,
                          dissolved_oxygen_downstream, pH_downstream,
                          turbidity_downstream, turbidity_downstream),
                .fns = ~is_na(.x)))

p_level_down <- data_downstream_imputed %>%
  select(Timestamp, level_downstream) %>%
  left_join(y = data_downstream_is.na %>%
              select(Timestamp, level_downstream),
            by = "Timestamp") %>%
  ggplot(aes(Timestamp, level_downstream.x, color = level_downstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_downstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, level_downstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"),
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("Level (m)")

p_temp_down <- data_downstream_imputed %>%
  select(Timestamp, temperature_downstream) %>%
  left_join(y = data_downstream_is.na %>%
              select(Timestamp, temperature_downstream),
            by = "Timestamp") %>%
  ggplot(aes(Timestamp, temperature_downstream.x,
             color = temperature_downstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_downstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, temperature_downstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"),
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab(expression('Temperature ('~degree*C*')'))

p_cond_down <- data_downstream_imputed %>%
  select(Timestamp, conductance_downstream) %>%
  left_join(y = data_downstream_is.na %>%
              select(Timestamp, conductance_downstream),
            by = "Timestamp") %>%
  ggplot(aes(Timestamp, conductance_downstream.x,
             color = conductance_downstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_downstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, conductance_downstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"),
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab(expression(Conductance~(mu*S*'/'*cm)))

p_turb_down <- data_downstream_imputed %>%
  select(Timestamp, turbidity_downstream) %>%
  left_join(y = data_downstream_is.na %>%
              select(Timestamp, turbidity_downstream),
            by = "Timestamp") %>%
  ggplot(aes(Timestamp, turbidity_downstream.x,
             color = turbidity_downstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_downstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, turbidity_downstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"),
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("turbidity (FNU)")

p_DO_down <- data_downstream_imputed %>%
  select(Timestamp, dissolved_oxygen_downstream) %>%
  left_join(y = data_downstream_is.na %>%
              select(Timestamp, dissolved_oxygen_downstream),
            by = "Timestamp") %>%
  ggplot(aes(Timestamp, dissolved_oxygen_downstream.x,
             color = dissolved_oxygen_downstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_downstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, dissolved_oxygen_downstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"),
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("DO (mg/l)")

p_pH_down <- data_downstream_imputed %>%
  select(Timestamp, pH_downstream) %>%
  left_join(y = data_downstream_is.na %>%
              select(Timestamp, pH_downstream),
            by = "Timestamp") %>%
  ggplot(aes(Timestamp, pH_downstream.x,
             color = pH_downstream.y)) +
  geom_point(size = 0.5) +
  geom_point(data = data_downstream_imputed %>%
               filter(Timestamp >= ymd("2019-10-24")
                      & Timestamp < ymd_hms("2019-10-25 12:00:00") |
                        Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00") |
                        Timestamp >= ymd_hms("2019-11-28 10:45:00")
                      & Timestamp < ymd_hms("2019-11-29 11:10:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, pH_downstream)) +
  scale_colour_manual(values = c("#000000", "#D55E00"),
                      breaks = c(0,1)) +
  theme(legend.position = "none") +
  ylab("pH")


ggpubr::ggarrange(p_turb_down, p_cond_down, p_DO_down, p_pH_down, p_level_down,
                  p_temp_down, ncol = 1)





###########################################################
###########################################################


## ---- lagTimeEstimation

# y_t: turbidity_downstream
# x_t: turbidity_upstream
# z_t: level_upstream, temperature_upstream

data_xz <- data_upstream_imputed %>% 
  select(Timestamp, level_upstream, temperature_upstream,
         turbidity_upstream)

data_y <- data_downstream_imputed %>% 
  select(Timestamp, turbidity_downstream)

data_lag_estim <- data_xz %>% 
  left_join(data_y, by = "Timestamp") %>% 
  mutate(turbidity_upstream = if_else(turbidity_upstream == 0, 
                                      as.numeric(NA), 
                                      turbidity_upstream))


## ---- Conditional-cross-correlation
mean <- c(5,5)
variance <- c(4,4)
correlation <- c(3,3)

# we choose level and temperature for lag time estimation in this analysis
cond_ccf <- data_lag_estim %>%
  conduits::conditional_ccf(x = turbidity_upstream,
                            y = turbidity_downstream,
                            z_numeric = c(level_upstream, 
                                          temperature_upstream),
                            k = 1:24,
                            knots_mean = list(x = mean, 
                                              y = mean),
                            knots_variance = list(x = variance, 
                                                  y = variance),
                            df_correlation = correlation)

## ---- Estimating dt
Estimate_dt <- cond_ccf %>% 
  estimate_dt(new_data = data_lag_estim, k_min = 1, 
              k_max = 24)

# 
# ## ----Figure (conditional-ccf)
# 
# # Comparing conditional ccf at lag dt and lag k for k=1,...,24
# 
# # 1. First compute conditionally normalised x_t and yt+dt conditional on zt
# # 2. Then compute x_t*yt+dt*
# # 3. Finally computing E[x_t*yt+dt*|z_t] by fitting a model - we choose 
# #    to fit the same model structure as we fitted for x_t*y_t+k
# 
# data <- Estimate_dt$data_dt
# 
# up <- c("turbidity_upstream",
#         "level_upstream", "temperature_upstream")
# down <- c("turbidity_downstream")
# lag <- "dt"
# 
# data <- data %>%
#   dplyr::select(.data$Timestamp, {{up}}, {{down}}, {{lag}})
# 
# L_min <- data %>% dplyr::pull({{lag}}) %>% min(na.rm = T)
# L_max <- data %>% dplyr::pull({{lag}}) %>% max(na.rm = T)
# 
# # computing lead timestamps at k=1,...,24
# Time <- data$Timestamp
# 
# Time_leads <- purrr::map(L_min:L_max, ~dplyr::lead(data$Timestamp, 
#                                                    n = .x))
# names(Time_leads) <- c(L_min:L_max)
# 
# Time_leads <- Time_leads %>%
#   tibble::as_tibble() %>%
#   dplyr::mutate(Timestamp = Time)
# 
# df_time_lead <- data %>%
#   dplyr::select(.data$Timestamp, {{lag}}) %>%
#   dplyr::left_join(Time_leads, by = "Timestamp")
# 
# df_ref_time_lead <- df_time_lead %>%
#   tidyr::pivot_longer(cols = -c(.data$Timestamp, {{lag}}),
#                       names_to = "lead", values_to = "time_d") %>%
#   dplyr::select(.data$Timestamp, .data$lead, .data$time_d) %>%
#   dplyr::rename(max_lag = .data$lead) %>%
#   dplyr::mutate(max_lag = as.double(.data$max_lag))
# 
# 
# 
# # Maping lags to get lead times corresponds to y_t+dt
# data_dt <- data %>%
#   dplyr::select(.data$Timestamp, {{up}}, {{lag}}) %>%
#   dplyr::rename(max_lag = {{lag}}) %>%
#   dplyr::inner_join(df_ref_time_lead, by = c("Timestamp", "max_lag"))
# 
# # joining downstream data with estimated time leads
# df_downstream <- data %>%
#   dplyr::select(.data$Timestamp, {{down}}) 
# 
# 
# data_dt <- data_dt %>% 
#   inner_join(df_downstream, c("time_d" = "Timestamp"))
# 
# data_dt <- tibble(Timestamp = Time) %>% 
#   full_join(data_dt) %>% 
#   rename("turbidity_downstream_dt" = turbidity_downstream)
# 
# 
# ##-- Computing conditional moments for x_t and y_t+dt
# # We use the same smoothing parameters as the ones we used before
# dim_mean_x <- cond_ccf$other$knots_mean$x
# dim_var_x <- cond_ccf$other$knots_var$x
# 
# dim_mean_y <- cond_ccf$other$knots_mean$y
# dim_var_y <- cond_ccf$other$knots_var$y
# 
# cond_moments_xt <- conditional_moments(data = data_dt, 
#                                        x = turbidity_upstream, 
#                                        z_numeric = c(level_upstream, 
#                                                      temperature_upstream),
#                                        knots_mean = dim_mean_x,
#                                        knots_variance = dim_var_x)
# 
# cond_moments_yt_dt <- conditional_moments(data = data_dt, 
#                                           x = turbidity_downstream_dt, 
#                                           z_numeric = c(level_upstream, 
#                                                         temperature_upstream),
#                                           knots_mean = dim_mean_y,
#                                           knots_variance = dim_var_y)
# 
# # conditional moment data for xt
# data_xt <- cond_moments_xt$data_conditional_moments %>% 
#   select(Timestamp, turbidity_upstream, E_turbidity_upstream,
#          Var_turbidity_upstream)
# 
# # conditional moment data for y_t+dt
# data_yt <- cond_moments_yt_dt$data_conditional_moments %>% 
#   select(Timestamp, turbidity_downstream_dt, E_turbidity_downstream_dt,
#          Var_turbidity_downstream_dt)
# 
# # combining all data and computing x_t*, y_t+dt* and xt*yt+dt*
# data_cond_moments <- data_dt %>% 
#   select(Timestamp, level_upstream, temperature_upstream) %>% 
#   list(data_xt, data_yt) %>% 
#   reduce(left_join, by = "Timestamp") %>% 
#   mutate(turbidity_upstream_star = (turbidity_upstream - E_turbidity_upstream)/sqrt(Var_turbidity_upstream),
#          turbidity_downstream_dt_star = (turbidity_downstream_dt - E_turbidity_downstream_dt)/sqrt(Var_turbidity_downstream_dt),
#          XY_star = turbidity_upstream_star*turbidity_downstream_dt_star)
# 
# # Now let's compute the conditional CCFs. Note that we use the same model 
# # that was fitted in cond_ccf_mod1
# 
# model_formula <- str_replace(cond_ccf$formula_gam,
#                              "XY",
#                              "XY_star") 
# df_correlation <- cond_ccf$other$df_correlation
# 
# # User defined link function
# corrlink <- function() {
#   ## link
#   linkfun <- function(mu) {log((1+mu)/(1-mu))}
#   ## inverse link
#   linkinv <- function(eta) {(exp(eta) - 1)/(exp(eta) + 1)}
#   ## derivative of invlink wrt eta
#   mu.eta <- function(eta) { 2*exp(eta)/(exp(eta) + 1)^2 }
#   valideta <- function(eta) TRUE
#   link <- "corrlink"
#   structure(list(linkfun = linkfun,
#                  linkinv = linkinv,
#                  mu.eta = mu.eta,
#                  valideta = valideta,
#                  name = link),
#             class = "link-glm")
# }
# 
# corrl <- corrlink()
# gam_fit_xydt <- stats::glm(formula = stats::as.formula(model_formula),
#                            data = data_cond_moments,
#                            family = stats::gaussian(link = corrl),
#                            start = rep(0,(sum(df_correlation))),
#                            control = stats::glm.control(maxit = 400))
# 
# 
# ccf_xydt <- stats::predict.glm(gam_fit_xydt, newdata = data_cond_moments, 
#                                type = "response")
# 
# data_cond_ccf <- cond_ccf$data_ccf %>% 
#   select(Timestamp, str_subset(names(cond_ccf$data_ccf),
#                                "k =")) %>% 
#   mutate("k = dt" = ccf_xydt) 
# 
# 
# # select(-turbidity_downstream, -turbidity_upstream,
# #        -level_upstream, -temperature_upstream,
# #        -conductance_upstream, -dissolved_oxygen_upstream)
# 
# # summary(data_cond_ccf)
# 
# 
# data_cond_ccf %>% 
#   pivot_longer(cols = -c(Timestamp,`k = dt`), names_to = "lag",
#                values_to = "ccf") %>% 
#   mutate(lag = factor(lag, levels = c(paste("k = ", L_min:L_max, 
#                                             sep = ""), "k = dt"))) %>% 
#   ggplot(aes(Timestamp, ccf, color = lag)) +
#   geom_line() +
#   geom_line(data = data_cond_ccf, aes(y = `k = dt`), color = "black") +
#   ylab("conditional ccf")
# 




## ---- Figure (vis-dt)


Plot_dt <- autoplot(object = Estimate_dt, m=1000, seed = 1989, 
                    interval = TRUE)
Plot_dt <- Plot_dt$plot_list
#--Adjusting the x axis in the dt plots

Plot_dt[[1]] <- Plot_dt[[1]] + 
  xlab(expression(Level~(m)))
Plot_dt[[2]] <- Plot_dt[[2]] + 
  xlab(expression('Temperature ('~degree*C*')'))

gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = Plot_dt, ncol = 2))


## ---- plot-dt-ccf-pval

insig_times <- Estimate_dt$data_dt %>% 
  filter(pvalue > 0.01) %>% 
  pull(Timestamp)

df1 <- Estimate_dt$data_dt %>% 
  select(Timestamp, dt, ccf, pvalue) %>% 
  pivot_longer(-Timestamp) 
  

df1 %>% 
  mutate(name = factor(name, levels = c("dt", "ccf", "pvalue"))) %>% 
  ggplot(aes(Timestamp, value)) +
  geom_line() +
  geom_point(data = df1 %>% 
              filter(Timestamp %in% insig_times), 
            aes(Timestamp, value), color = "red", size = 0.3) +
  facet_wrap(~name, ncol = 1, scales = "free_y")


#############################################
##-- Replacing insignificant dt --##
############################################

# Finding the lag that gives maximum cross correlation 
# between between xt and y_t+k for k = 1,...,24

CCF <- forecast::Ccf(x = data_xy$turbidity_upstream, 
                     y = data_xy$turbidity_downstream, 
                     lag.max = 24, plot=F)

df_ccf <- tibble(lag = tail(CCF$lag, 24), ccf = tail(CCF$acf, 24))

max_lag_ccf <- df_ccf %>% 
  dplyr::slice(which.max(ccf)) %>% 
  pull(lag) %>% 
  as.numeric()

data_dt <- Estimate_dt$data_dt %>% 
  mutate(dt = if_else(pvalue > 0.01, max_lag_ccf, dt))




############################################################
############################################################

## ----Computing upstream variables with estimated lags

data_dt <- data_dt %>% 
  select(Timestamp, turbidity_downstream, dt, ccf) 
data_upstream_new <- list(data_upstream_imputed, data_dt) %>% 
  reduce(left_join, by = "Timestamp")

data_upstream_new_dt <- data_upstream_new %>% 
  map_lags(up = c(level_upstream, conductance_upstream,
                  dissolved_oxygen_upstream,
                  turbidity_upstream, temperature_upstream,
                  pH_upstream), 
           down = c(turbidity_downstream), 
           lag = dt)


data_lagged_upstream <- data_upstream_new_dt %>% 
  select(-turbidity_downstream)

save(data_lagged_upstream, file = here::here("weird_river_data_paper", "data", "data_lagged_upstream.rda"))
save(data_downstream_imputed, file = here::here("weird_river_data_paper", "data", "data_downstream_imputed.rda"))






