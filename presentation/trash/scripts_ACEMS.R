library(tidyverse)
library(lubridate)
library(mgcv)
library(patchwork)

## ---- tsplots

load(here::here("data", "PRIN_5min_flagged.rda"))
load(here::here("data", "PRIN_5min_cleaned.rda"))

PRIN_5min_flagged <-  PRIN_5min_flagged %>% 
  mutate(Timestamp = ymd_hms(roundedTimestamp)) %>% 
  filter(Timestamp >= ymd("2019-10-01") &
           Timestamp < ymd("2020-01-01"))

PRIN_5min_cleaned <-  PRIN_5min_cleaned %>% 
  mutate(Timestamp = ymd_hms(roundedTimestamp)) %>% 
  filter(Timestamp >= ymd("2019-10-01") &
           Timestamp < ymd("2020-01-01")) %>% 
  rename("level" = surfacewaterElevMean,
         "temperature" = surfWaterTempMean,
         "conductance" = specificConductance,
         "dissolved_oxygen" = dissolvedOxygen)

PRIN_data <- PRIN_5min_flagged %>% 
  select(Timestamp, site, surfacewaterElevMean, surfWaterTempMean, 
         specificConductance, turbidity, dissolvedOxygen) %>% 
  rename("level" = surfacewaterElevMean,
         "temperature" = surfWaterTempMean,
         "conductance" = specificConductance,
         "dissolved_oxygen" = dissolvedOxygen)


p_turb <- PRIN_data %>% 
  select(Timestamp, turbidity, site) %>% 
  ggplot(aes(Timestamp, turbidity, color = site)) +
  geom_line() +
  scale_color_manual(values = c("#999999", "#E69F00"),
                     breaks = c("down", "up")) +
  ylab("Turbidity (FNU)")

p_cond <- PRIN_data %>% 
  select(Timestamp, conductance, site) %>% 
  ggplot(aes(Timestamp, conductance, color = site)) +
  geom_line() +
  scale_color_manual(values = c("#999999", "#E69F00"),
                     breaks = c("down", "up")) +
  ylab(expression(Conductance~(mu*S*'/'*cm)))

p_do <- PRIN_data %>% 
  select(Timestamp, dissolved_oxygen, site) %>% 
  ggplot(aes(Timestamp, dissolved_oxygen, color = site)) +
  geom_line() +
  scale_color_manual(values = c("#999999", "#E69F00"),
                     breaks = c("down", "up")) +
  ylab("DO (mg/l)")

p_level <- PRIN_data %>% 
  select(Timestamp, level, site) %>% 
  ggplot(aes(Timestamp, level, color = site)) +
  geom_line() +
  scale_color_manual(values = c("#999999", "#E69F00"),
                     breaks = c("down", "up")) +
  ylab("Level (m)")

p_temp <- PRIN_data %>% 
  select(Timestamp, temperature, site) %>% 
  ggplot(aes(Timestamp, temperature, color = site)) +
  geom_line() +
  scale_color_manual(values = c("#999999", "#E69F00"),
                     breaks = c("down", "up")) +
  ylab(expression('Temperature ('~degree*C*')')) +
  theme(legend.position = "bottom")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p_temp)

grid.arrange(arrangeGrob(p_cond + theme(legend.position = "none"),
                         p_do + theme(legend.position = "none"), 
                         p_level + theme(legend.position = "none"),
                         p_temp + theme(legend.position = "none"),
                         p_turb + theme(legend.position = "none"),
                         ncol = 2), 
             mylegend, heights=c(10, 1))










## ---- turbplot

Time_turb_down_anomaly <- PRIN_5min_flagged %>% 
  filter(site == "down",
         turbidityAnomalyFlag == 1) %>% 
  pull(Timestamp)

plot_turb_down <- PRIN_data %>%
  filter(site == "down") %>% 
  select(Timestamp, turbidity) %>% 
  ggplot(aes(Timestamp, turbidity)) +
  geom_point(alpha = 0.3, size = 0.3) +
  geom_line() +
  geom_point(data = PRIN_data %>% filter(site == "down",
                                         Timestamp %in% Time_turb_down_anomaly),
             aes(Timestamp, turbidity), color = "red", size = 1) +
  geom_point(data = PRIN_data %>%
               filter(Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, turbidity)) +
  ylab("Turbidity (FNU)") +
  ggtitle("downstream")


plot_turb_up <- PRIN_5min_cleaned %>%
  filter(site == "up") %>% 
  select(Timestamp, turbidity) %>% 
  ggplot(aes(Timestamp, turbidity)) +
  geom_point(alpha = 0.3, size = 0.3) +
  geom_line() +
  geom_point(data = PRIN_5min_cleaned %>%
               filter(Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, turbidity)) +
  ylab("Turbidity (FNU)")+
  ggtitle("upstream")


plot_level <- PRIN_5min_cleaned %>%
  filter(site == "up") %>% 
  select(Timestamp, level) %>% 
  ggplot(aes(Timestamp, level)) +
  geom_point(alpha = 0.3, size = 0.3) +
  geom_line() +
  geom_point(data = PRIN_5min_cleaned %>%
               filter(site == "up",
                      Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, level)) +
  ylab("Level (m)") +
  ggtitle("upstream")

  
plot_cond <- PRIN_5min_cleaned %>%
  filter(site == "up") %>% 
  select(Timestamp, conductance) %>% 
  ggplot(aes(Timestamp, conductance)) +
  geom_point(alpha = 0.3, size = 0.3) +
  geom_line() +
  geom_point(data = PRIN_5min_cleaned %>%
               filter(site == "up",
                      Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, conductance)) +
  ylab(expression(Conductance~(mu*S*'/'*cm))) +
  ggtitle("upstream")

plot_temp <- PRIN_5min_cleaned %>%
  filter(site == "up") %>% 
  select(Timestamp, temperature) %>% 
  ggplot(aes(Timestamp, temperature)) +
  geom_point(alpha = 0.3, size = 0.3) +
  geom_line() +
  geom_point(data = PRIN_5min_cleaned %>%
               filter(site == "up",
                      Timestamp >= ymd("2019-11-07")
                      & Timestamp < ymd_hms("2019-11-09 00:00:00")),
             color = "#009E73", size = 0.5,
             aes(Timestamp, temperature)) +
  ylab(expression('Temperature ('~degree*C*')'))+
  ggtitle("upstream")


grid.arrange(plot_turb_down + theme(axis.title.x = element_blank()), 
             plot_turb_up + theme(axis.title.x = element_blank()), 
             plot_cond + theme(axis.title.x = element_blank()),
             plot_level + theme(axis.title.x = element_blank()), 
             plot_temp, 
             ncol = 1)





