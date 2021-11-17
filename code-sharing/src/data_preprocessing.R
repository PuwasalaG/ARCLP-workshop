# See weird_river_data_paper/data-raw/Identifying_QF.Rmd for identifying anomalies
# in variables

## ---- load-waq-upANDdown

load(here::here("weird_river_data_paper", "data-raw", "flagged-data",
                "waq_down.rda"))
load(here::here("weird_river_data_paper", "data-raw", "flagged-data",
                "waq_up.rda"))

# View(waq_down)
# names(waq_down)

waq_up <- waq_up %>%
  filter(Timestamp >= ymd("2019-10-01") 
         & Timestamp < ymd("2020-01-01")) %>% 
  mutate(turbidityAnomalyFinalFlag = if_else(turbidityAnomalyFlag == 1 
                                             | turbidity_wiperAnomalyFlag == 1,
                                             1, turbidityAnomalyFlag)) 
waq_down <- waq_down %>%
  filter(Timestamp >= ymd("2019-10-01") 
         & Timestamp < ymd("2020-01-01")) %>% 
  mutate(turbidityAnomalyFinalFlag = if_else(turbidityAnomalyFlag == 1 
                                             | turbidity_wiperAnomalyFlag == 1,
                                             1, turbidityAnomalyFlag)) 


## ---- rawData-timeplots-upstream
p_turb_up <- waq_up %>% 
  select(Timestamp, turbidity, turbidityAnomalyFlag) %>% 
  ggplot(aes(Timestamp, turbidity)) +
  geom_point(size = 0.3) +
  geom_line() +
  geom_point(data = waq_up %>% 
               filter(turbidityAnomalyFlag == 1),
             aes(Timestamp, turbidity, color = "red")) +
  theme(legend.position = "none") +
  ylab("Turbidity (FNU)")

p_cond_up <- waq_up %>% 
  select(Timestamp, specificConductance, conductanceAnomalyFlag) %>% 
  ggplot(aes(Timestamp, specificConductance)) +
  geom_point(size = 0.3) +
  geom_line() +
  ylab(expression(Conductance~(mu*S*'/'*cm)))
  # geom_point(data = waq_up %>% 
  #              filter(conductanceAnomalyFlag == 1),
  #            aes(Timestamp, specificConductance, color = "red"))

p_do_up <- waq_up %>% 
  select(Timestamp, dissolvedOxygen, dissolvedOxygenAnomalyFlag) %>% 
  ggplot(aes(Timestamp, dissolvedOxygen)) +
  geom_point(size = 0.3) +
  geom_line() +
  ylab("DO (mg/l)")
  # geom_point(data = waq_up %>% 
  #              filter(dissolvedOxygenAnomalyFlag == 1),
  #            aes(Timestamp, dissolvedOxygen, color = "red"))

p_pH_up <- waq_up %>% 
  select(Timestamp, pH, pHAnomalyFlag) %>% 
  ggplot(aes(Timestamp, pH)) +
  geom_line() +
  geom_point(size = 0.3) +
  ylab("pH")
  # geom_point(data = waq_up %>% 
  #              filter(pHAnomalyFlag == 1),
  #            aes(Timestamp, pH, color = "red"))


ggpubr::ggarrange(p_turb_up, p_cond_up, p_do_up, p_pH_up, ncol = 1)


## ---- rawData-timeplots-downstream

p_turb_down <- waq_down %>% 
  select(Timestamp, turbidity, turbidityAnomalyFlag) %>% 
  ggplot(aes(Timestamp, turbidity)) +
  geom_point(size = 0.3) +
  geom_line() +
  geom_point(data = waq_down %>% 
               filter(turbidityAnomalyFlag == 1),
             aes(Timestamp, turbidity, color = "red")) +
  theme(legend.position = "none") +
  ylab("Turbidity (FNU)")

p_cond_down <- waq_down %>% 
  select(Timestamp, specificConductance, conductanceAnomalyFlag) %>% 
  ggplot(aes(Timestamp, specificConductance)) +
  geom_point(size = 0.3) +
  geom_line() +
  theme(legend.position = "none") +
  ylab(expression(Conductance~(mu*S*'/'*cm))) +
  geom_point(data = waq_down %>%
             filter(conductanceAnomalyFlag == 1),
             aes(Timestamp, specificConductance, color = "red"))

p_do_down <- waq_down %>% 
  select(Timestamp, dissolvedOxygen, dissolvedOxygenAnomalyFlag) %>% 
  ggplot(aes(Timestamp, dissolvedOxygen)) +
  geom_point(size = 0.3) +
  geom_line() +
  ylab("DO (mg/l)")
# geom_point(data = waq_down %>% 
#              filter(dissolvedOxygenAnomalyFlag == 1),
#            aes(Timestamp, dissolvedOxygen, color = "red"))

p_pH_down <- waq_down %>% 
  select(Timestamp, pH, pHAnomalyFlag) %>% 
  ggplot(aes(Timestamp, pH)) +
  geom_point(size = 0.3) +
  geom_line() +
  ylab("pH")
# geom_point(data = waq_down %>% 
#              filter(pHAnomalyFlag == 1),
#            aes(Timestamp, pH, color = "red"))


ggpubr::ggarrange(p_turb_down, p_cond_down, p_do_down, p_pH_down, 
                  ncol = 1)



## ---- wiper-anomaly

waq_down %>% 
  filter(Timestamp >= ymd("2019-10-14") & Timestamp < ymd("2019-10-15")) %>% 
  mutate(wiper_anomaly_flag = factor(turbidity_wiperAnomalyFlag,
                                     levels = c(0,1), 
                                     labels = c("valid data", "wiper anomaly"))) %>% 
  ggplot(aes(Timestamp, turbidity, 
             color = wiper_anomaly_flag)) +
  geom_point() +
  theme(legend.position = "none") +
  ylab("Turbidity (FNU)") +
  scale_color_manual(values = c("black", "red"))

  