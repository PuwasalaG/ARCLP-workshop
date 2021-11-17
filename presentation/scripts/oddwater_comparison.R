## ----setup, include=FALSE------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  dev = "png",
  cache = TRUE,
  warning = FALSE,
  message = FALSE)
library(tidyverse)
library(lubridate)
library(fpp3)
library(patchwork)
library(conduits)



## ----loadoutPRINData-----------------------------
load(here::here("Analysis", "oddwater-comparison", "data",
                "PRIN_5min_flagged.rda"))

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

p1<- PRIN_5min_turb %>%
  ggplot(aes(x= Timestamp, y= turbidity )) + 
  geom_line() +
  geom_point(aes(colour = flag)) +
  scale_color_manual(values = c("o" = "black", "1" = "red"))

print(p1)



## ----oct-----------------------------------------

turb_oct <- PRIN_5min_turb %>%
   filter(Timestamp >= ymd("2019-10-01") & 
           Timestamp < ymd("2020-01-01"))

p1<- turb_oct %>%
  ggplot(aes(x= Timestamp, y= turbidity )) + 
  geom_line() +
  geom_point(aes(colour = flag)) +
  scale_color_manual(values = c("o" = "black", "1" = "red"))

print(p1)


## ------------------------------------------------
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

plot_trans <- function(data, col = "turbidity")
{
  p<- data %>%
  ggplot(aes(x= Timestamp, y= data[[col]] )) + 
  geom_line() +
  geom_point(aes(colour = flag)) +
  scale_color_manual(values = c("o" = "black", "1" = "red"))+
  ylab(col)

return(p)

}





## ----oddwateroct---------------------------------
data <- turb_oct %>%
  select(-flag) %>%
  tidyr::drop_na() 
trans_data <-transform_data1(data)
full_data <- left_join(trans_data, turb_oct %>% select(-turbidity), by = "Timestamp")


p1 <- plot_trans(data = full_data, col = "turbidity")
print(p1)
p2 <- plot_trans(data = full_data, col = "log_turbidity")
print(p2)
p3<- plot_trans(data = full_data, col = "diff_turbidity")
print(p3)
p4<- plot_trans(data = full_data, col = "difflog_log_turbidity")
print(p4)
p5<- plot_trans(data = full_data, col = "rc_turbidity")
print(p5)
p6<- plot_trans(data = full_data, col = "ratio_turbidity")
print(p6)
p7<- plot_trans(data = full_data, col = "rdifflog_turbidity")
print(p7)
p8<- plot_trans(data = full_data, col = "rdiff_turbidity")
print(p8)


(p1|p2)/
(p3|p4)/
(p5|p6)/
(p7|p8)


## ------------------------------------------------
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
out_rc <- stray::find_HDoutliers(data[["rc_turbidity"]], p= p, alpha = alpha, k=k)

turb_data <- data %>%
  select(Timestamp, turbidity, flag) %>%
  mutate(type = out_rc$type)

code_trub <- mapply(
    code, turb_data[["flag"]],
    turb_data[["type"]]
  )  

turb_data <- turb_data %>%
  mutate(code  = code_trub)


sum_tub <- list(p=p, alpha = alpha, k=k, FN = sum(turb_data$code=="FN") ,  FP = sum(turb_data$code=="FP") ,  TN = sum(turb_data$code=="TN"),  TP = sum(turb_data$code=="TP"))

return(sum_tub %>% rbind() %>% as_tibble())


}

sumtab <- find_performance(data =full_data, p=0.5, alpha=0.01 ,k=1)

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.5, alpha=0.01 ,k=5))

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.75, alpha=0.01 ,k=1))

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.75, alpha=0.01 ,k=5))

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.75, alpha=0.05 ,k=1))
sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.75, alpha=0.05 ,k=5))

sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.80, alpha=0.05 ,k=5))
sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.80, alpha=0.05 ,k=10))
sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.80, alpha=0.1 ,k=5))
sumtab <-bind_rows(sumtab, find_performance(data =full_data, p=0.80, alpha=0.1 ,k=10))

knitr::kable(sumtab)



## ------------------------------------------------
out_rc <- stray::find_HDoutliers(full_data$rc_turbidity, p= 0.75, alpha = 0.05, k=5)

turb_data <- full_data %>%
  select(Timestamp, turbidity, flag) %>%
  mutate(type = out_rc$type)

code_trub <- mapply(
    code, turb_data[["flag"]],
    turb_data[["type"]]
  )  

table(code_trub)
turb_data <- turb_data %>%
  mutate(code  = code_trub)

p<- turb_data %>%
  ggplot(aes(x= Timestamp, y= turbidity )) + 
  geom_line() +
  geom_point(aes(colour = code)) 
  #scale_color_manual(values = c("o" = "black", "1" = "red"))+

print(p)

plotly::ggplotly(p)


