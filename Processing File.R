## CO2 d13C data from Picarro 
## 03092018 

setwd("C:/Users/Chris Wilson/Desktop/BRU Respiration Logs/Picarro_Logs/03092018")


library(lubridate); library(ggplot2); library(dplyr)

CO2_log <- read.table("03092018_p1c1.dat",header = T)

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


CO2_log_plateaus <- CO2_log %>% filter(CO2_deriv < 0.5 & CO2_deriv > - 0.5, X12CO2 > 200)



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
setwd("C:/Users/Chris Wilson/Desktop/BRU Respiration Logs/Field_Logs/03152018")
list.files()
field_CO2 <- read.table("p5_a1.txt", header = TRUE)
str(field_CO2);

field_CO2$sec <- seq(1,1530,1);
 

?apply

filenames2 <- list.files(pattern = "*.txt")
flux_df_list <- lapply(filenames2, read.table, header = T)
# OK, I have a list, but it is unnamed. 
names(flux_df_list) <- substr(filenames2,1,5)

flux_est <- rep(0,length(names(flux_df_list)));
flux_se <- rep(0,length(flux_est)); 

for(i in 1:length(flux_est)){

field_CO2 <- flux_df_list[[i]];
field_CO2$sec <- seq(1,length(field_CO2$Time.H.M.S.),1);
field_CO2$CO2_smooth <- smooth.spline(field_CO2$sec, field_CO2$CO2.ppm.)[[2]]

plot_trace <- ggplot(field_CO2, aes(x=sec, y = CO2_smooth)) + geom_line()
 print(plot_trace)
 
ambient <- 400
flux_df <- field_CO2 %>% filter(CO2.ppm. > ambient - 15, CO2.ppm. < ambient + 15, 
                                sec > 100) # will need to figure out better time filter 

flux_deriv <- diff(flux_df$CO2.ppm.)/diff(flux_df$sec)

plot(flux_df$sec[-1],flux_deriv)
flux_df$flux_deriv <- c(0,flux_deriv)


flux_df2 <- flux_df %>% filter(flux_deriv > -2.5, flux_deriv < 2.5)

measure <- rep(0, length(flux_df2$flux_deriv));
measure[1] <- 1;

for(i in 2:length(measure)){
  if(flux_df2$sec[i] <= flux_df2$sec[i-1]+20){
    measure[i] <- measure[i-1];
  } else{
    measure[i] <- measure[i-1] + 1;
  }
  if(measure[i]>3) warning('measure[i] > 3!', i)
}

# BADOOM! 

flux_df2$measure <- measure

library(broom); library(dplyr); library(ggplot2)

ggplot(flux_df2, aes(x=sec,y=flux_deriv, color = measure)) + geom_point()

# summarizing differentials
fluxes_df <- flux_df2 %>% group_by(measure) %>% summarize(m_flux = mean(flux_deriv,na.rm = T),
                                                                     se_flux  = sd(flux_deriv,na.rm = T))
# fitting lines 
fluxes_lm <- flux_df2 %>% group_by(measure) %>% do(fitSlope = lm(CO2.ppm. ~ sec, data = .))
lm_coefs <- tidy(fluxes_lm, fitSlope)
print(lm_coefs)

plot_fit <- ggplot() + geom_point(data = flux_df2, aes(x = sec, y = CO2.ppm.,color = measure)) +
  geom_abline(aes(intercept = lm_coefs$estimate[1], slope = lm_coefs$estimate[2]), color = "red") +
  geom_abline(aes(intercept = lm_coefs$estimate[3], slope = lm_coefs$estimate[4]), color = "red") +
  geom_abline(aes(intercept = lm_coefs$estimate[5], slope = lm_coefs$estimate[6]), color = "red") 
  
print(plot_fit)

flux_est[i] <- coef(summary(flux_model))[2,1]
flux_se[i] <- coef(summary(flux_model))[2,2]

}



### Playing around with smoothers 
spl <- smooth.spline(CO2_log$t[-893], y=CO2_log$Delta_Raw_iCO2[-893])
pred <- predict(spl)

plot (CO2_log$t, CO2_log$Delta_Raw_iCO2, log="xy")
lines(pred, col=2)

ycs.prime <- diff(CO2_log$X12CO2[-893])/diff(CO2_log$t[-893])
pred.prime <- predict(spl, deriv=1)


C_pore <- 450 
C_atm <- rep(0,3*10^2)
C_atm[1] <- 350
k1 <- 0.03
k2 <- 1

for(i in 2:length(C_atm)){
  C_atm[i] <- C_atm[i-1] + k1*(C_pore - C_atm[i-1]) + k2;
}

plot(C_atm)
time <- seq(1,length(C_atm),1)
which(C_atm > 450)

c_diff <- diff(C_atm)/diff(time)

plot(c_diff)




y <- rnorm(5,0,1)

param_proposal <- seq(-5,5,0.1)
log_like <- rep(0,length(param_proposal))
for(i in 1:length(param_proposal)){
log_like[i] <- -sum(log(dnorm(y,param_proposal[i],1)))
}


plot(param_proposal,exp(-log_like))
sum(exp(-log_like))



