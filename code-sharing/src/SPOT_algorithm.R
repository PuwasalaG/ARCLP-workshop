# Functions for SPOT algorithm

find_threshold_POT <- function(X, threshold = NULL, 
                               threshold_prob = 0.9,
                               level = 0.001, 
                               cmax = FALSE, r = 1){
  
  X <- na.omit(X) %>% abs() %>% as.numeric()
  
  # calculate threshold
  t1 <- threshold
  if(is.null(t1)){
    t1 <- quantile(X, probs = min(threshold_prob,(1-level)))
  }
  
  n <- length(X)
  Nt <- sum(X>t1)
  
  # estimate scale and shape parameters of GPD using POT
  fit_gpd <- evd::fpot(X, threshold = t1, std.err = F, 
                       cmax = cmax, r = r)
  
  # compute outlier threshold
  scale_gpd <- fit_gpd$estimate[1]
  shape_gpd <- fit_gpd$estimate[2]
  thresh <- evd::qgpd(p = 1-level, loc = t1, scale = scale_gpd,
                      shape = shape_gpd)
  
  
  k <- sum(X>thresh)
  
  return(list(k = k, outlier_threshold = thresh, 
              initial_threshold = t1))
  
}


# Streaming Peak-Over-Threshold

find_outliers_SPOT <- function(residuals,
                               n,
                               initial_threshold = NULL,
                               initial_threshold_prob = 0.9,
                               level = 0.001, cmax = FALSE, r = 1){
  
  # residuals to calibrate outlier threshold
  X_calib <- residuals %>% slice(1:n) %>% drop_na()
  X_calib_resid <- X_calib$residuals
  
  # streaming resiudals
  X_stream <- residuals %>% slice(n+1:n()) %>% drop_na()
  
  
  # Computing outlier threshold from POT
  POT <- find_threshold_POT(X = X_calib_resid,
                            threshold_prob = initial_threshold_prob,
                            level = level, cmax = cmax, r = r)
  z_q <- POT$outlier_threshold
  t1 <- POT$initial_threshold
  
  # identifying outliers in the calibration set
  
  calib_ol_timestamp <- X_calib %>% 
    filter(residuals > z_q|residuals < -z_q) %>% 
    pull(Timestamp)
  
  # removing outliers in the calibration set at re-calibrating z_q
  if(!is.null(calib_ol_timestamp)){
    X_calib <- X_calib %>%
      mutate(residuals = if_else(Timestamp %in% calib_ol_timestamp,
                                 as.numeric(NA), residuals))
    POT <- find_threshold_POT(X = X_calib$residuals,
                              threshold_prob = initial_threshold_prob,
                              level = level, cmax = cmax, r = r)
    z_q <- POT$outlier_threshold
    X_calib_resid <- X_calib$residuals %>% na.omit()
    
    calib_ol_timestamp2 <- X_calib %>% 
      filter(residuals > z_q|residuals < -z_q) %>% 
      pull(Timestamp)
    
    calib_ol_timestamp <- calib_ol_timestamp %>% 
      append(calib_ol_timestamp2)
  }
  
  k <- length(X_calib_resid)
  Nt <- sum(X_calib_resid>t1)
  list_ol_timestamp <- numeric()
  
  for (i in 1:nrow(X_stream)) {
    
    x_i <- X_stream %>% magrittr::extract(i,2) %>% as.numeric()
    
    if(x_i > z_q|x_i < -z_q){
      list_ol_timestamp[i] <- X_stream %>% 
        magrittr::extract(i,1) %>% 
        as.numeric()
    } else if(x_i > t1|x_i < -t1){
      X_calib_resid <- c(X_calib_resid, x_i)
      X1 <- X_calib_resid %>% abs()
      Nt <- Nt+1
      k <- k+1
      
      # estimate scale and shape parameters of GPD using POT
      fit_gpd <- evd::fpot(X1, threshold = t1, std.err = F,
                           cmax = cmax, r = r)
      
      # compute outlier threshold
      scale_gpd <- fit_gpd$estimate[1]
      shape_gpd <- fit_gpd$estimate[2]
      z_q <- evd::qgpd(p = 1-level, loc = t1, scale = scale_gpd,
                       shape = shape_gpd)
      
    } else{
      X_calib_resid <- c(X_calib_resid, x_i)
    }
    
    
  }
  
  list_ol_timestamp <- list_ol_timestamp %>% 
    na.omit() %>% 
    append(calib_ol_timestamp) %>% 
    lubridate::as_datetime() %>% 
    sort()
  
  residuals <- residuals %>% 
    mutate(AnomalyFlag_SPOT = if_else(Timestamp %in% list_ol_timestamp,
                                            "outlier", "typical"),
           AnomalyFlag_SPOT = if_else(is.na(residuals), 
                                            as.character(NA),
                                      AnomalyFlag_SPOT))
  
  return(residuals)
}

## ---- function-to-compute-confusion-matrix

compute_confusion <- function(data, residuals, name, n, q, t_prob,
                              cmax = FALSE, r = 1){
  resids <- data %>%
    select(Timestamp, {{residuals}}) %>%
    rename("residuals" = {{residuals}}) %>%
    mutate(residuals = as.numeric(residuals))
  
  resids <- find_outliers_SPOT(residuals = resids, n = n,
                               initial_threshold_prob = t_prob,
                               level = q, cmax = cmax, r = r)
  
  resids <- data %>%
    select(Timestamp, turbidity_downstream, turbidityAnomalyFlag) %>%
    left_join(resids)
  
  resids <- resids %>%
    mutate(confusion_matrix = case_when(turbidityAnomalyFlag
                                        == "typical"
                                        & AnomalyFlag_SPOT ==
                                          "outlier" ~ "FP",
                                        turbidityAnomalyFlag
                                        == "outlier"
                                        & AnomalyFlag_SPOT ==
                                          "typical" ~ "FN",
                                        turbidityAnomalyFlag
                                        == "outlier"
                                        & AnomalyFlag_SPOT ==
                                          "outlier" ~ "TP",
                                        turbidityAnomalyFlag
                                        == "typical"
                                        & AnomalyFlag_SPOT ==
                                          "typical" ~ "TN",
                                        turbidityAnomalyFlag
                                        == "outlier"
                                        & is.na(AnomalyFlag_SPOT)
                                        ~ "FN",
                                        TRUE ~
                                          as.character(NA)),
           confusion_matrix = factor(confusion_matrix,
                                     levels = c("TP", "TN", "FP", "FN",
                                                "NA"))) %>%
    select(Timestamp, AnomalyFlag_SPOT, confusion_matrix)
  
  
  return(resids)
}

