## CO2 d13C data from Picarro 
## 03092018 ---- 

setwd("C:/Users/Chris Wilson/Desktop/BRU Respiration Logs/Picarro_Logs/03152018")
library(lubridate); library(ggplot2); library(dplyr)

library(plyr)

filenames <- list.files(pattern = "*.dat")
ldf <- lapply(filenames, read.table, header = T)
# OK, I have a list, but it is unnamed. 
names(ldf) <- substr(filenames,1,13)
str(ldf[[1]])


# Split element 4 into element 4 and 5 
ldf_4 <- ldf[[4]][which(CO2_log$t < 950),]
ldf_5 <- ldf[[4]][which(CO2_log$t >= 950),]
ldf[[4]] <- ldf_4
ldf2 <- ldf
ldf2 <- append(ldf2, list(ldf_5), 4)


names(ldf2)[5] <- "03092018_p1e1"

# Can use the ddply functions here...eg...
lapply(ldf, ddply, .(country), summarise, spend  = sum(spend),
       trials = sum(trials))


keeling_function(ldf, filenames)




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


















# Creating vector to store delta point estimates
delta_13C_Reco <- rep(0,length(filenames));
delta_13C_Reco_SE <- rep(0, length(filenames)); 
r_squared_vals <- rep(0,length(filenames)); 

for(i in 1:(length(filenames))){
CO2_log <- ldf[[i]]

CO2_log$time <- hms(CO2_log$TIME)
CO2_log$t <- time_length(CO2_log$time) - min(time_length(CO2_log$time))
CO2_log$CO2_smooth <- smooth.spline(CO2_log$t, CO2_log$X12CO2)[[2]]



ggplot(CO2_log, aes(x = t)) + ylab("CO2 (ppm)") +
  theme_bw() + geom_point(aes(y = CO2_smooth)) + 
  geom_point(aes(x = t, y = -50*Delta_Raw_iCO2), color = "red") +
  scale_y_continuous(limits = c(0,1000))

ggplot(CO2_log, aes(x = t)) + ylab("CO2 (ppm)") +
  theme_bw() + geom_point(aes(y = Delta_Raw_iCO2)) + scale_y_continuous(limits = c(-20,20))



# Evaluating derivative criterion for identifying plateaus 
CO2_deriv <- diff(CO2_log$CO2_smooth)/diff(CO2_log$t)

## Optional plotting
deriv_df <- data.frame(x = CO2_log$t[-1], y = CO2_deriv)
ggplot(deriv_df, aes(x=x,y=y)) + geom_point() + scale_y_continuous(limits = c(-1,1))

## Suggests that tolerances should be |dCO2/dt| around 0.25  
CO2_log$CO2_deriv <- c(0,CO2_deriv); 
CO2_log_plateaus <- CO2_log %>% filter(CO2_deriv < 0.3 & CO2_deriv > - 0.3, X12CO2 > 150)

## Optional plotting
plottyness <- ggplot(CO2_log_plateaus, aes(x = t)) + ylab("CO2 (ppm)") +
  theme_bw() + geom_point(aes(y = CO2_smooth)) + 
  geom_point(aes(x = t, y = -50*Delta_Raw_iCO2), color = "red") +
  scale_y_continuous(limits = c(0,1500));
print(plottyness); 

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


keeling_data <- CO2_log_plateaus %>% group_by(syringe) %>% summarize(mCO2 = mean(X12CO2,na.rm = T),
                                                                     md13C = mean(Delta_Raw_iCO2,na.rm = T)) %>%
  mutate(inv_CO2 = 1/mCO2);

plot(keeling_data$inv_CO2, keeling_data$md13C)
collar_Reco_d13C <- lm(md13C ~ inv_CO2, keeling_data);
summary(collar_Reco_d13C)

delta_13C_Reco[i] <- coef(summary(collar_Reco_d13C))[1,1]
delta_13C_Reco_SE[i] <- coef(summary(collar_Reco_d13C))[1,2]
r_squared_vals[i] <- summary(collar_Reco_d13C)$r.squared

print(i);
}











CO2_log <- read.table("03092018_p1c1.dat",header = T)





CO2_log <- ldf2[[4]]

CO2_log$time <- hms(CO2_log$TIME)
CO2_log$t <- time_length(CO2_log$time) - min(time_length(CO2_log$time))
CO2_log$CO2_smooth <- smooth.spline(CO2_log$t, CO2_log$X12CO2)[[2]]

ggplot(CO2_log, aes(x = t)) + ylab("CO2 (ppm)") +
  theme_bw() + geom_point(aes(y = CO2_smooth)) + 
  geom_point(aes(x = t, y = -50*Delta_Raw_iCO2), color = "red") +
  scale_y_continuous(limits = c(0,1000))


ggplot(CO2_log, aes(x = t, y = Delta_Raw_iCO2)) + geom_point() +
  geom_line(color = "black",size = 0.5) + scale_y_continuous(limits = c(-20,5)) +
  ylab("d13C") + theme_bw() + scale_x_continuous(breaks = seq(0,1900,100));


ggplot(CO2_log, aes(x = t, y = Delta_30s_iCO2)) + geom_point() +
  geom_line(color = "black",size = 0.5) + scale_y_continuous(limits = c(-20,5)) +
  ylab("d13C") + theme_bw() + scale_x_continuous(breaks = seq(0,1900,100));


# Evaluating derivative criterion for identifying plateaus 
CO2_deriv <- diff(CO2_log$CO2_smooth)/diff(CO2_log$t)
plot(CO2_log$t[-1],CO2_deriv)
deriv_df <- data.frame(x = CO2_log$t[-1], y = CO2_deriv)
ggplot(deriv_df, aes(x=x,y=y)) + geom_point() + scale_y_continuous(limits = c(-1,1))
## Suggests that tolerances should be |dCO2/dt| around 0.25  
CO2_log$CO2_deriv <- c(0,CO2_deriv); 

CO2_log_plateaus <- CO2_log %>% filter(CO2_deriv < 0.3 & CO2_deriv > - 0.3, X12CO2 > 200)



ggplot(CO2_log_plateaus, aes(x = t)) + ylab("CO2 (ppm)") +
  theme_bw() + geom_point(aes(y = CO2_smooth)) + 
  geom_point(aes(x = t, y = -50*Delta_Raw_iCO2), color = "red") +
  scale_y_continuous(limits = c(0,1000))

## BADOOM! 
# Now need to separate syringes 

syringe <- rep(0, length(CO2_log_plateaus$t));
syringe[1] <- 1
for(i in 2:length(syringe)){
  if(CO2_log_plateaus$t[i] <= CO2_log_plateaus$t[i-1]+30){
    syringe[i] <- syringe[i-1];
  } else{
    syringe[i] <- syringe[i-1] + 1;
  }
  if(syringe[i]>4) warning('Syringe[i] > 4!')
}

# BADOOM! 

CO2_log_plateaus$syringe <- syringe 


keeling_data <- CO2_log_plateaus %>% group_by(syringe) %>% summarize(mCO2 = mean(X12CO2,na.rm = T),
                                                                     md13C = mean(Delta_Raw_iCO2,na.rm = T)) %>%
  mutate(inv_CO2 = 1/mCO2);


plot(keeling_data$inv_CO2, keeling_data$md13C)
collar_Reco_d13C <- lm(md13C ~ inv_CO2, keeling_data);
summary(collar_Reco_d13C)













# Intervals 125-250, 375-480, 625-710, 875-950 

library(dplyr)

CO2_log %>% filter(t > 800 & t < 900) %>% summarise(meand13C = mean(Delta_Raw_iCO2),
                                                    meanCO2 = mean(X12CO2))
### E.g. manual extraction of Rhizoma peanut 
peanut_delta_p2 <- c(-8.935, -8.2126, -7.22, -7.114);
peanut_CO2_p2 <- c(354,396,469,492);

peanut_CO2_p2_inv <- 1/peanut_CO2_p2;

peanut_source_delta_p2 <- lm(peanut_delta_p2 ~ peanut_CO2_p2_inv)
summary(peanut_source_delta_p2); # Intercept is -2.25 [0.2475]

# Background delta 
delta_back <- function(delta_a, c_a, delta_s = -2.25, c_b = 354){
  delta_b = (delta_a - delta_s)*(c_a/c_b) + delta_s;
  return(delta_b);
}
delta_back(peanut_delta_p2, peanut_CO2_p2); # -8.9 


peanut_delta_p1 <- c(-2.286, -1.10, -0.21, 0.27);
peanut_CO2_p1 <- c(423,442,496,519);

peanut_CO2_p1_inv <- 1/peanut_CO2_p2;

peanut_source_delta_p1 <- lm(peanut_delta_p1 ~ peanut_CO2_p1_inv)
summary(peanut_source_delta_p1); # Intercept is 6.41 [0.65]

delta_back(peanut_delta_p1, peanut_CO2_p1, 6.41, 423); # -1.54


# Field Logs 

setwd("C:/Users/Chris Wilson/Desktop/BRU Respiration Logs/Field_Logs/03092018")
list.files()
field_CO2 <- read.table("p1_b1.txt", header = TRUE)
str(field_CO2);
field_CO2$sec <- seq(1,length(field_CO2$Time.H.M.S.),1);


### Playing around with smoothers 
spl <- smooth.spline(CO2_log$t[-893], y=CO2_log$Delta_Raw_iCO2[-893])
pred <- predict(spl)

plot (CO2_log$t, CO2_log$Delta_Raw_iCO2, log="xy")
lines(pred, col=2)

ycs.prime <- diff(CO2_log$X12CO2[-893])/diff(CO2_log$t[-893])
pred.prime <- predict(spl, deriv=1)


16*24/60
