


keeling_function <- function(ldf, filenames){
  
  # Creating vector to store delta point estimates
  delta_13C_Reco <- rep(0,length(filenames));
  delta_13C_Reco_SE <- rep(0, length(filenames)); 
  r_squared_vals <- rep(0,length(filenames)); 
  
  #pdf(paste(substr(filenames[i],1,13),".pdf"))
  pdf("Master_pdf.pdf")
  
  #(length(filenames))
  for(i in 1:(length(filenames))){
    
    CO2_log <- ldf[[i]]
    
    # Computing time and CO2 smoothed 
    CO2_log$time <- hms(CO2_log$TIME)
    CO2_log$t <- time_length(CO2_log$time) - min(time_length(CO2_log$time))
    CO2_log$CO2_smooth <- smooth.spline(CO2_log$t, CO2_log$X12CO2)[[2]]
    
    ggplot(CO2_log, aes(x=t, y = Delta_Raw_iCO2)) + geom_point()
    
    # Evaluating derivative criterion for identifying plateaus 
    CO2_deriv <- diff(CO2_log$CO2_smooth)/diff(CO2_log$t)
    
    ## Optional plotting
    deriv_df <- data.frame(x = CO2_log$t[-1], y = CO2_deriv)
    #  ggplot(deriv_df, aes(x=x,y=y)) + geom_point() + scale_y_continuous(limits = c(-1,1))
    
    ## Suggests that tolerances should be |dCO2/dt| around 0.25  
    CO2_log$CO2_deriv <- c(0,CO2_deriv); 
    CO2_log_plateaus <- CO2_log %>% filter(CO2_deriv < 0.3 & CO2_deriv > - 0.3, X12CO2 > 150)
    
    syringe <- rep(0, length(CO2_log_plateaus$t));
    syringe[1] <- 1
    for(j in 2:length(syringe)){
      if(CO2_log_plateaus$t[j] <= CO2_log_plateaus$t[j-1]+30){
        syringe[j] <- syringe[j-1];
      } else{
        syringe[j] <- syringe[j-1] + 1;
      }
      if(syringe[j]>4) warning('Syringe[j] > 4!',i)
    }
    
    # BADOOM! 
    CO2_log_plateaus$syringe <- syringe 
    
    CO2_log_plateaus_filter <- CO2_log_plateaus %>% group_by(syringe) %>%
      filter(quantile(Delta_Raw_iCO2, 0.9) > Delta_Raw_iCO2) %>%
      filter(quantile(Delta_Raw_iCO2, 0.1) < Delta_Raw_iCO2)
    
    
    
    ## Optional plotting
    plottyness <- ggplot(CO2_log_plateaus_filter, aes(x = t)) + ylab("CO2 (ppm)") +
      theme_bw() + geom_point(aes(y = CO2_smooth)) + 
      geom_point(aes(x = t, y = -50*Delta_Raw_iCO2), color = "red") +
      scale_y_continuous(limits = c(0,1500));
    
    print(plottyness); 
    
    
    keeling_data <- CO2_log_plateaus_filter %>% group_by(syringe) %>% summarize(mCO2 = mean(X12CO2,na.rm = T),
                                                                                md13C = mean(Delta_Raw_iCO2,na.rm = T)) %>% mutate(inv_CO2 = 1/mCO2);
    
    
    ## pdf(paste(substr(filenames[1],1,13),".pdf"));
    plot(keeling_data$inv_CO2, keeling_data$md13C);
    
    
    collar_Reco_d13C <- lm(md13C ~ inv_CO2, keeling_data);
    ## summary(collar_Reco_d13C)
    
    delta_13C_Reco[i] <- coef(summary(collar_Reco_d13C))[1,1]
    delta_13C_Reco_SE[i] <- coef(summary(collar_Reco_d13C))[1,2]
    r_squared_vals[i] <- summary(collar_Reco_d13C)$r.squared
    
    print(i);
  }
  dev.off();
  return(list(delta_13C = delta_13C_Reco, delta_13C_SE = delta_13C_Reco_SE, 
              rSquared = r_squared_vals))
  
}
