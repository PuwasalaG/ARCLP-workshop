### Supplemental Materials (Data and scripts)
### Anomaly Detection in River Networks

## ---- load
library(tidyverse)
library(lubridate)
library(fpp3)
library(patchwork)
library(conduits)


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

load(here::here("weird_river_data_paper", "data", "PRIN_5min_flagged.rda"))
PRIN_5min_flagged <- PRIN_5min_flagged %>%
  rename(
    Timestamp = roundedTimestamp,
    level = surfacewaterElevMean,
    conductance = specificConductance,
    dissolved_oxygen = dissolvedOxygen,
    temperature = surfWaterTempMean
  ) %>%
  mutate(Timestamp = ymd_hms(Timestamp))

## ---- tsanomalyplotPRINData
# Pivot data values from wide to long
data_val <- PRIN_5min_flagged %>%
  select(Timestamp, conductance, dissolved_oxygen,
       pH, turbidity, chlorophyll, fDOM,
       level, temperature, site ) %>%
  as_tsibble(index = Timestamp, key = site) %>%
  pivot_longer(cols = conductance:temperature) 

# Pivot data flags from wide to long
data_flag <- PRIN_5min_flagged %>%
  select(Timestamp, conductanceAnomalyFlag, 
         dissolvedOxygenAnomalyFlag,
         pHAnomalyFlag, turbidityAnomalyFlag,
         chlorophyllAnomalyFlag, fDOMAnomalyFlag,
         site ) %>%
  rename(
    conductance = conductanceAnomalyFlag,
    dissolved_oxygen = dissolvedOxygenAnomalyFlag,
    pH = pHAnomalyFlag,
    turbidity = turbidityAnomalyFlag,
    chlorophyll = chlorophyllAnomalyFlag,
    fDOM = fDOMAnomalyFlag
  ) %>%
  as_tsibble(index = Timestamp, key = site) %>%
  pivot_longer(cols = conductance:fDOM)

# Combine values and flag data

data<- left_join(data_val, data_flag, by = 
                   c("Timestamp", "site", "name" ) ) %>%
      rename(value = value.x,
             flag = value.y)

data_outliers <- data %>%
  filter(flag == 1)
  

p <- data %>%
  ggplot(aes(x = Timestamp, y = value, color = site)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("down" = "black", "up" = "blue")) +
  theme(legend.position='bottom') +
  geom_point(data = data_outliers, aes(x = Timestamp, y = value), color = "red" )

print(p)

## ---- tsplotUPstream

data_outliers <- data %>%
  filter(flag == 1, site == "up")

p1 <- data %>%
  filter(site == "up") %>%
  ggplot(aes(x = Timestamp, y = value)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  theme(legend.position='bottom') +
  geom_point(data = data_outliers, aes(x = Timestamp, y = value), color = "red" )

print(p1)

## ---- tsplotDownstream

data_outliers <- data %>%
  filter(flag == 1, site == "down")

p2 <- data %>%
  filter(site == "down") %>%
  ggplot(aes(x = Timestamp, y = value)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  theme(legend.position='bottom') +
  geom_point(data = data_outliers, aes(x = Timestamp, y = value), color = "red" )

print(p2)

## ---- tsplotPRINData

PRIN_5min_new <- PRIN_5min_cleaned %>%
  filter(Timestamp >= ymd("2019-10-01")
  & Timestamp < ymd("2020-01-01")) %>%
  select(str_subset(names(PRIN_5min_cleaned), "Anomaly", negate = T)) # Check this step

PRIN_5min_new <- PRIN_5min_new %>%
  as_tsibble(index = Timestamp, key = site)

p <- PRIN_5min_new %>%
  pivot_longer(cols = conductance:temperature) %>%
  ggplot(aes(x = Timestamp, y = value, color = site)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("down" = "black", "up" = "red")) +
  theme(legend.position='bottom')
  
print(p)

## ---- filterDataForAnalysis

# Data prepreration
# Spreading the values for up and down sites

data_z <- PRIN_5min_new %>%
  filter(site == "up") %>%
  select(-site, -fDOM, -turbidity, -chlorophyll, -pH) %>%
  rename(
    "level_upstream" = level,
    "temperature_upstream" = temperature,
    "conductance_upstream" = conductance,
    "dissolved_oxygen_upstream" = dissolved_oxygen
  )


data_xy <- PRIN_5min_new %>%
  select(Timestamp, site, turbidity) %>%
  spread(key = site, value = turbidity) %>%
  rename(
    "turbidity_downstream" = down,
    "turbidity_upstream" = up
  )

data_PRIN_sub <- data_xy %>%
  left_join(data_z, by = "Timestamp") %>%
  mutate(turbidity_upstream = if_else(turbidity_upstream == 0,
    as.numeric(NA),
    turbidity_upstream
  ))

## ---- pairwiseScatterplotsUp

p1 <- data_PRIN_sub %>%
  select(-turbidity_downstream) %>%
  pivot_longer(
    cols = conductance_upstream:temperature_upstream,
    names_to = "upstream_var", values_to = "value"
  ) %>%
  ggplot(aes(value, turbidity_upstream)) +
  geom_point() +
  facet_wrap(~upstream_var, scales = "free_x") +
  theme(aspect.ratio = 1)
print(p1)


## ---- pairwiseScatterplotsDown
p2 <- data_PRIN_sub %>%
  select(-turbidity_upstream) %>%
  pivot_longer(
    cols = conductance_upstream:temperature_upstream,
    names_to = "upstream_var", values_to = "value"
  ) %>%
  ggplot(aes(value, turbidity_downstream)) +
  geom_point() +
  facet_wrap(~upstream_var, scales = "free_x") +
  theme(aspect.ratio = 1)

print(p2)

## ---- scatterpair
p <- data_PRIN_sub %>%
  select(
    turbidity_downstream, turbidity_upstream, conductance_upstream,
    level_upstream, dissolved_oxygen_upstream,
    temperature_upstream
  ) %>%
  GGally::ggpairs()
print(p)

## ---- CCF
# cross-correlation functions

p1 <- data_PRIN_sub %>%
  CCF(turbidity_downstream, turbidity_upstream) %>%
  autoplot() +
  ggtitle("Turbidity downstream and turbidity upstream")

p2 <- data_PRIN_sub %>%
  CCF(turbidity_downstream, conductance_upstream) %>%
  autoplot() +
  ggtitle("Turbidity downstream and conductance upstream")

p3 <- data_PRIN_sub %>%
  CCF(turbidity_downstream, dissolved_oxygen_upstream) %>%
  autoplot() +
  ggtitle("Turbidity downstream and  dissolved oxygen upstream")

p4 <- data_PRIN_sub %>%
  CCF(turbidity_downstream, level_upstream) %>%
  autoplot() +
  ggtitle("Turbidity downstream and level upstream")

p5 <- data_PRIN_sub %>%
  CCF(turbidity_downstream, temperature_upstream) %>%
  autoplot() +
  ggtitle("Turbidity downstream and temperature upstream")

(p1 | p2) /
  (p3 | p4) /
  p5

## ---- VAR
var1_fit <- MTS::VARMA(data, p = 1, q = 0, include.mean = FALSE, details = F)
print(var1_fit)
MTS::ccm(var1_fit$residuals, lags = 25)

## ---- cond_ccf

cond_ccf_mod1 <- data_PRIN_sub %>%
  as_tibble() %>%
  conduits::conditional_ccf(
    x = turbidity_upstream,
    y = turbidity_downstream,
    z_numeric = c(
      level_upstream,
      conductance_upstream,
      temperature_upstream
    ),
    k = 1:24,
    knots_mean = list(
      x = c(8, 8, 8),
      y = c(8, 8, 8)
    ),
    knots_variance = list(
      x = c(7, 7, 7),
      y = c(7, 7, 7)
    ),
    df_correlation = c(4, 2, 4)
  )

p1 <- cond_ccf_mod1 %>% autoplot()
print(p1)

# Estimate_dt_mod1 <- cond_ccf_mod1 %>%
#  estimate_dt(new_data = data_PRIN_sub, k_min = 1,
#              k_max = 30)

Estimate_dt_mod1 <- cond_ccf_mod1 %>%
  estimate_dt(
    new_data = data_PRIN_sub, k_min = 1,
    k_max = 24
  )

Estimate_dt_mod1$data_dt %>% summary()

p2 <- Estimate_dt_mod1 %>% autoplot()
print(p2)

## ---- compute_upstream_with_dt

# data_dt_mod1 <- Estimate_dt_mod1$data_dt %>%
#  select(Timestamp, dt, ccf) %>%
#  rename("dt_mod1" = dt,
#        "ccf_mod1" = ccf)
# data_PRIN_new <- left_join(data_PRIN_sub, data_dt_mod1, by = "Timestamp")
# data_PRIN_new <- list(data_PRIN_sub, data_dt_mod1) %>%
#   reduce(left_join, by = "Timestamp")

data_PRIN_new_dt_mod1 <- Estimate_dt_mod1$data_dt %>%
  map_lags(
    up = c(
      level_upstream, conductance_upstream,
      dissolved_oxygen_upstream,
      turbidity_upstream, temperature_upstream
    ),
    down = c(turbidity_downstream),
    lag = dt
  )

# https://github.com/tidyverse/tibble/blob/master/R/subsetting.R
# Error: Can't use NA as column index in a tibble for assignment.
head(data_PRIN_new_dt_mod1)

## ---- Modeling_downstream_turbidity

load(here::here("weird_river_data_paper", "data", "PRIN_5min_flagged.rda"))
PRIN_5min_flagged <- PRIN_5min_flagged %>%
  rename(
    Timestamp = roundedTimestamp,
    level = surfacewaterElevMean,
    conductance = specificConductance,
    dissolved_oxygen = dissolvedOxygen,
    temperature = surfWaterTempMean
  ) %>%
  mutate(Timestamp = ymd_hms(Timestamp))

PRIN_5min_flagged_turbidity <- PRIN_5min_flagged %>%
  filter(
    Timestamp >= ymd("2019-10-01") &
      Timestamp < ymd("2020-01-01"),
    site == "down"
  ) %>%
  select(Timestamp, turbidity, turbidityAnomalyFlag) %>%
  mutate(turbidity = if_else(turbidity <= 0, as.numeric(NA),
    turbidity
  )) %>%
  rename("turbidity_downstream" = turbidity)

data_mod_turb <- data_PRIN_new_dt_mod1 %>%
  select(-turbidity_downstream) %>%
  left_join(PRIN_5min_flagged_turbidity)
