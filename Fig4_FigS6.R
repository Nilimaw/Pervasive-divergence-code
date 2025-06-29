##Rcode for Fig4_CD
##First model from Greenfield 2006 supplementary paper for AHA1
##Second model from Greenfield 2006 supplementary paper for GUK1
##Fit for normalized mean residue molar ellipticity
library(dplyr)
library(tidyr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library(ggtext)
library(showtext)
library(nls2)
library(investr)
library(confintr)
library(stats)
library(nlstools)

theme_set(theme_classic())
setwd("C:/Users/Nilima/Desktop/Fig4-20250318T220527Z-001/Fig4")
##Read Tm data 
#S. cer
C1 <- read.csv("Tm_CGUK1_R1_edited.csv", header = T)
C2 <- read.csv("Tm_CGUK1_R2_edited.csv", header = T)
C3 <- read.csv("Tm_CGUK1_R3_edited.csv", header = T)

#S. uv
U1 <- read.csv("Tm_UGUK1_R1_edited.csv", header = T)
U2 <- read.csv("Tm_UGUK1_R2_edited.csv", header = T)
U3 <- read.csv("Tm_UGUK1_R3_edited.csv", header = T)

C1 <- C1 %>% drop_na()
C1$species <- "cer"
C1$Rep <- 1
C2 <- C2 %>% drop_na()
C2$species <- "cer"
C2$Rep <- 2
colnames(C2)[2] <- "Theta"
colnames(C2)[1] <- "Temperature"
C3 <- C3 %>% drop_na()
C3$species <- "cer"
C3$Rep <- 3
colnames(C3)[2] <- "Theta"

U1 <- U1 %>% drop_na()
U1$species <- "uv"
U1$Rep <- 1
U2 <- U2 %>% drop_na()
U2$species <- "uv"
U2$Rep <- 2
U3 <- U3 %>% drop_na()
U3$species <- "uv"
U3$Rep <- 3
colnames(U3)[2] <- "Theta"

##Calculate molar ellipticity and mean residue molar ellipticity
##https://www.photophysics.com/faqs/methods-techniques/cd-units-faqs/
#https://www.bioline.com/media/calculator/01_04.html
#Molar concentration of S. cer
#conc. by Bradford assay of sample = 0.25mg/mL
#Mol. weight of S. cer GUK1 = 20637.37 Da
#Mol. weight of S. uv GUK1 = 20560.33 Da
#Length of GUK1 = 187
#Path length of cuvette = 1mm
x = 0.25/20637.37
cG_M = 0.00001211395 # (Molar)
y = 0.25/20560.33
uG_M = 0.00001215934
C1$molar_ellip <- 100*C1$Theta/(0.1*cG_M*1000)
C1$mean_residue_molar_ellip <- C1$molar_ellip/187
C2$molar_ellip <- 100*C2$Theta/(0.1*cG_M*1000)
C2$mean_residue_molar_ellip <- C2$molar_ellip/187
C3$molar_ellip <- 100*C3$Theta/(0.1*cG_M*1000)
C3$mean_residue_molar_ellip <- C3$molar_ellip/187

U1$molar_ellip <- 100*U1$Theta/(0.1*uG_M*1000)
U1$mean_residue_molar_ellip <- U1$molar_ellip/187
U2$molar_ellip <- 100*U2$Theta/(0.1*uG_M*1000)
U2$mean_residue_molar_ellip <- U2$molar_ellip/187
U3$molar_ellip <- 100*U3$Theta/(0.1*uG_M*1000)
U3$mean_residue_molar_ellip <- U3$molar_ellip/187


##Take average of three replicates for GUK1 for Temperature, Theta, molar ellipticity and mean residual molar ellipticity
cdf <- cbind(C1,C2,C3)
colnames(cdf)
colnames(cdf)[8] <- "Temp.2"
colnames(cdf)[15] <- "Temp.3"
colnames(cdf)[9] <- "Theta.2"
colnames(cdf)[16] <- "Theta.3"
colnames(cdf)[13] <- "molar_ellip_2"
colnames(cdf)[20] <- "molar_ellip_3"
colnames(cdf)[14] <- "mean_residue_molar_ellip_2"
colnames(cdf)[21] <- "mean_residue_molar_ellip_3"

cdf <- cdf %>%
  select(Temperature,Theta,molar_ellip,mean_residue_molar_ellip, Temp.2,Theta.2,molar_ellip_2,mean_residue_molar_ellip_2, Temp.3,Theta.3, molar_ellip_3, mean_residue_molar_ellip_3)
cdf <- cdf %>%
  mutate(cer_Temp_avg = (Temperature + Temp.2 + Temp.3)/3) %>%
  mutate(cer_Theta_avg = (Theta + Theta.2 + Theta.3)/3) %>%
  mutate(cer_molar_ellip_avergae = (molar_ellip + molar_ellip_2 + molar_ellip_3)/3) %>%
  mutate(cer_mean_residue_avg = (mean_residue_molar_ellip + mean_residue_molar_ellip_2 + mean_residue_molar_ellip_3)/3)


udf <- cbind(U1,U2,U3)
colnames(udf)
colnames(udf)[8] <- "Temp.2"
colnames(udf)[15] <- "Temp.3"
colnames(udf)[9] <- "Theta.2"
colnames(udf)[16] <- "Theta.3"
colnames(udf)[13] <- "molar_ellip_2"
colnames(udf)[20] <- "molar_ellip_3"
colnames(udf)[14] <- "mean_residue_molar_ellip_2"
colnames(udf)[21] <- "mean_residue_molar_ellip_3"

udf <- udf %>%
  select(Temperature,Theta,molar_ellip,mean_residue_molar_ellip, Temp.2,Theta.2,molar_ellip_2,mean_residue_molar_ellip_2, Temp.3,Theta.3, molar_ellip_3, mean_residue_molar_ellip_3)

udf <- udf %>%
  mutate(uv_Temp_avg = (Temperature + Temp.2 + Temp.3)/3) %>%
  mutate(uv_Theta_avg = (Theta + Theta.2 + Theta.3)/3) %>%
  mutate(uv_molar_ellip_avergae = (molar_ellip + molar_ellip_2 + molar_ellip_3)/3) %>%
  mutate(uv_mean_residue_avg = (mean_residue_molar_ellip + mean_residue_molar_ellip_2 + mean_residue_molar_ellip_3)/3)

cdf <- cdf %>%
  select(cer_Temp_avg,cer_Theta_avg,cer_molar_ellip_avergae,cer_mean_residue_avg)
cdf$species <- "cer"
colnames(cdf) <- c("Temp", "Theta", "molar_ellip_average", "mean_residue_avg", "species")

udf <- udf %>%
  select(uv_Temp_avg,uv_Theta_avg, uv_molar_ellip_avergae,uv_mean_residue_avg)
udf$species <- "uv"
colnames(udf) <- c("Temp", "Theta", "molar_ellip_average", "mean_residue_avg", "species")

##Normalize mean residue molar ellipticity by dividing my min(mean residue moalr ellipticity)
no <- min(cdf$mean_residue_avg)
cdf$norm_mrme <- cdf$mean_residue_avg/no

no <- min(udf$mean_residue_avg)
udf$norm_mrme <- udf$mean_residue_avg/no


####Fit with Greenfield equation
##Cer_GUK1 - Full

##Implementation of the second set of equations:
#II. A two-state transition of a monomer between folded and unfolded forms with correcting the 
#data for pre- and post-transition linear changes in ellipticity as a function of temperature.
##Parameters
h=-20 #starting enthalpy 
m=314.15/1000 #starting TM in Kelvin
u= min(cdf$norm_mrme) #Mean Residue Ellipticity 100% folded 
l= max(cdf$norm_mrme) #Mean Residue Ellipticity 100% unfolded 
u1=100 #linear correction folded function of temp 
l1=-45 #linear correction unfolded function of temp
paras <- list(m = 0.31415, h=-20, u = 0.2598773, l = 1, u1 = 100, l1 = -45)


#Equation
#((exp((h/(1.987*Temp))*((Temp/m)-1)))/(1+exp((h/(1.987*Temp))*((Temp/m)-1))))*((u+(u1*Temp))-(l+(l1*Temp)))+ (l+(l1*Temp)) 

#Scale temperature data to facilitate curve fitting - convert celsius to kelvin and then divide by 1000
tcdf <- cdf %>%
  select(Temp,norm_mrme)
tcdf$Temp <- tcdf$Temp+273.15
tcdf$Temp <- tcdf$Temp/1000

##Model fit with nls
model_tcdf <- nls(norm_mrme ~ ((exp((h/(1.987*Temp))*((Temp/m)-1)))/(1+exp((h/(1.987*Temp))*((Temp/m)-1))))*((u+(u1*Temp))-(l+(l1*Temp)))+ (l+(l1*Temp)), start = paras, data = tcdf, control = nls.control(maxiter = 100))
model_tcdf
##Estimate CIs for parametrs 
confint(model_tcdf, 'm', level = 0.95)



##Uv_GUK1
##Implementation of the first set of equations:
##I. A two-state transition of a monomer from a folded to unfolded state. This treatment assumes 
##that the heat capacity of the folded and unfolded states are equal.
##Parameters
h=-200 #starting enthalpy in cal/mol 
m=333.15/1000 #starting TM in Kelvin. 
l = max(udf$norm_mrme) #mean residue ellipticity of unfolded protein 
u = min(udf$norm_mrme) #mean residue ellipticity of 100% folded helical protein 
paras <- list(m = 0.33315, h=-200, u =0.2045367, l =1, u1 = 53, l1 = 16)

###Scale temperature data to facilitate curve fitting - convert celsius to kelvin and then divide by 1000

tudf <- udf %>%
  select(Temp,norm_mrme)
tudf$Temp <- tudf$Temp+273.15
tudf$Temp <- tudf$Temp/1000

##Model
model_tudf <- nls(norm_mrme ~ ((exp((h/(1.987*Temp))*((Temp/m)-1)))/(1+exp((h/(1.987*Temp))*((Temp/m)-1))))*((u+(u1*Temp))-(l+(l1*Temp)))+ (l+(l1*Temp)), start = paras, data = tudf, control = nls.control(maxiter = 100))
model_tudf

confint2(model_tudf, level = 0.95)


##Predicted values from model
cdf_pre <- data.frame(predFit(model_tcdf, interval = "confidence", level= 0.95))
cdf_pre_og <- data.frame(cbind(tcdf,cdf_pre))
cdf_pre_og$Temp <- cdf_pre_og$Temp*1000-273.15
cdf_pre_og$species <- "cer"

udf_pre <- data.frame(predFit(model_tudf, interval = "confidence", level= 0.95))
udf_pre_og <- data.frame(cbind(tudf,udf_pre))
udf_pre_og$Temp <- udf_pre_og$Temp*1000 - 273.15
udf_pre_og$species <- "uv"


##Plot observed values and fitted values
ndf <- rbind(cdf_pre_og,udf_pre_og)
p1 <- ggplot( data = ndf, aes(x=Temp, y=norm_mrme, color=species)) +
  geom_point(aes(color=species), size = 1)+
  #geom_point(shape = 1,size = 1, color = "grey")+
  geom_line(aes(x = Temp, y=fit))+
  geom_ribbon(data=ndf, aes(x=Temp, ymin=lwr, ymax=upr, fill=species), alpha=0.5, inherit.aes=F)+
  geom_vline(xintercept = 51.05, color = "#DA6C42FF",alpha=0.5)+
  geom_vline(xintercept = 45.45, color = "#225BB2FF",alpha=0.5)+
  scale_color_manual(labels = c("S. cerevisiae", "S. uvarum"),values=c('#DA6C42FF','#2166ac'))+
  ylab("Normalized secondary structure")+
  xlab("Temperature")+
  ggtitle("GUK1") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title = element_blank())
p1





##Wavelength scans at different temperature for S. cer GUK1
wv_cer <- read.csv("CGUK1_rep2_wv_edited.csv", header = T)
colnames(wv_cer) <- c("wavelength","4C","20C","30C","40C","50C","60C","75C")

wv_cer <- wv_cer %>%
  pivot_longer(!wavelength, names_to = "Temp", values_to = "Theta")
wv_cer$wavelength <- as.numeric(wv_cer$wavelength)
wv_cer$Theta <- as.numeric(wv_cer$Theta)
wv_cer$Temp <- factor(wv_cer$Temp, levels = c("4C",  "20C", "30C", "40C", "50C", "60C", "75C"))

#Calculate mean residue molar ellipticity
wv_cer$molar_ellip <- 100*wv_cer$Theta/(0.1*cG_M*1000)
wv_cer$mean_residue_molar_ellip <- wv_cer$molar_ellip/187

c1 <- ggplot( data = wv_cer, aes(x=wavelength, y=mean_residue_molar_ellip, group = Temp, color = Temp)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(n = 7, name = "Reds"),labels = c(paste("4","\u00b0","C", sep = ""),paste("20","\u00b0","C", sep = ""),paste("30","\u00b0","C", sep = ""),paste("40","\u00b0","C", sep = ""),paste("50","\u00b0","C", sep = ""),paste("60","\u00b0","C", sep = ""),paste("75","\u00b0","C", sep = "")))+
  ylab("Mean residue molar ellipticity")+
  xlab("Wavelength (nm)")+
  ggtitle("GUK1") + 
  ylim(-9000,4000)+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
c1


##Wavelength scans at different temperature for S. uv GUK1
wv_uv <- read.csv("uguk1_rep1_wv.csv", header = T)
colnames(wv_uv) <- c("wavelength","4C","20C","30C","40C","50C","60C","75C")

wv_uv <- wv_uv %>%
  pivot_longer(!wavelength, names_to = "Temp", values_to = "Theta")
wv_uv$wavelength <- as.numeric(wv_uv$wavelength)
wv_uv$Theta <- as.numeric(wv_uv$Theta)

wv_uv$Temp <- factor(wv_uv$Temp, levels = c("4C",  "20C", "30C", "40C", "50C", "60C", "75C"))

##Calculate mean residue molar ellipticity
wv_uv$molar_ellip <- 100*wv_uv$Theta/(0.1*uG_M*1000)
wv_uv$mean_residue_molar_ellip <- wv_uv$molar_ellip/187


u1 <- ggplot( data = wv_uv, aes(x=wavelength, y=mean_residue_molar_ellip, group = Temp, color = Temp)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(n = 7, name = "Blues"),labels = c(paste("4","\u00b0","C", sep = ""),paste("20","\u00b0","C", sep = ""),paste("30","\u00b0","C", sep = ""),paste("40","\u00b0","C", sep = ""),paste("50","\u00b0","C", sep = ""),paste("60","\u00b0","C", sep = ""),paste("75","\u00b0","C", sep = "")))+
  ylab("Mean residue molar ellipticity")+
  xlab("Wavelength (nm)")+
  ggtitle("GUK1") + 
  ylim(-9000,4000)+
  #theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))

u1


##AHA1 CD data analysis
CA1 <- read.csv("cAHA1_rep_tm_edited.csv", header = T)
UA1 <- read.csv("uaHA1_rep1_tm_edited.csv", header = T)
CA1 <- CA1 %>% drop_na()
CA1$species <- "cer"
CA1$Rep <- 1

UA1 <- UA1 %>% drop_na()
UA1$species <- "uv"
UA1$Rep <- 1

#Molar concentration of S. cer
#conc. by BCA assay of sample:
#S. cer AHA1 = 0.719mg/ml
#S. uv AHA1 = 1.289mg/ml
#Mol. weight of S. cer AHA1 = 39435.53
#Mol. weight of S. uv AHA1 = 39089.05
#Length of S. cer AHA1 = 350
#Length of S. uv AHA1 = 347
#Path length of cuvette = 1mm
x = 0.719/39435.53
cA_M = 0.00001823229

y= 1.289/39089.05
uA_M = 0.00003297599
CA1$molar_ellip <- 100*CA1$Theta/(0.1*cA_M*1000)
CA1$mean_residue_molar_ellip <- CA1$molar_ellip/350
UA1$molar_ellip <- 100*UA1$Theta/(0.1*uA_M*1000)
UA1$mean_residue_molar_ellip <- UA1$molar_ellip/347



##Uv_AHA1
##Implementation of the first set of equations:
##I. A two-state transition of a monomer from a folded to unfolded state. This treatment assumes 
##that the heat capacity of the folded and unfolded states are equal.

#Normalize mean residue molar ellipticity by its min value
no <- min(UA1$mean_residue_molar_ellip)
UA1$norm_mrme <- UA1$mean_residue_molar_ellip/no

###Scale temperature data to facilitate curve fitting - convert celsius to kelvin and then divide by 1000
tUA <- UA1 %>%
  select(Temperature,norm_mrme)
tUA$Temperature <- (tUA$Temperature + 273.15)/1000

##Parameters
h=-200 #starting enthalpy in cal/mol 
tm= 323.15/1000 #starting TM in Kelvin 
l = max(UA1$norm_mrme)  #mean residue ellipticity of unfolded protein 
u = min(UA1$norm_mrme) #mean residue ellipticity of 100% folded helical protein 



##Model
paras <- list(m = 0.32315, h=-200, u = 0.1449745, l = 1 )
model_ua1 <- nls(norm_mrme ~ ((u-l)*((exp((h/(1.987*Temperature))*((Temperature/m)-1)))/(1+(exp((h/(1.987*Temperature))*((Temperature/m)-1))))))+ (l) , start = paras, data = tUA)
model_ua1

confint2(model_ua1,level = 0.95)



##CAHA1
##Implementation of the first set of equations:
##I. A two-state transition of a monomer from a folded to unfolded state. This treatment assumes 
##that the heat capacity of the folded and unfolded states are equal.

###Normalize mean residue molar ellipticity by its min value
no <- min(CA1$mean_residue_molar_ellip)
CA1$norm_mrme <- CA1$mean_residue_molar_ellip/no

#Scale temperature data to facilitate curve fitting - convert celsius to kelvin and then divide by 1000
tCA <- CA1 %>%
  select(Temperature, norm_mrme)
tCA$Temperature <- (tCA$Temperature + 273.15)/1000

##Parameters
h=-200 #starting enthalpy in cal/mol 
tm=333.15/1000 #starting TM in Kelvin. 
l = max(CA1$norm_mrme) #mean residue ellipticity of unfolded protein
u = min(CA1$norm_mrme) #mean residue ellipticity of 100% folded helical protein 



#Model
paras <- list(m = 0.33315, h=-200, u =0.274847, l =1)
model_ca1 <- nls(norm_mrme ~ ((u-l)*((exp((h/(1.987*Temperature))*((Temperature/m)-1)))/(1+(exp((h/(1.987*Temperature))*((Temperature/m)-1))))))+ (l) , start = paras,data = tCA)
model_ca1
confint2(model_ca1,level = 0.95)


#Predict values 
cadf_pre <- data.frame(predFit(model_ca1, interval = "confidence", level= 0.95))
cadf_pre_og <- data.frame(cbind(tCA,cadf_pre))
cadf_pre_og$Temperature <- cadf_pre_og$Temp*1000-273.15
cadf_pre_og$species <- "cer"

uadf_pre <- data.frame(predFit(model_ua1, interval = "confidence", level= 0.95))
uadf_pre_og <- data.frame(cbind(tUA,uadf_pre))
uadf_pre_og$Temperature <- uadf_pre_og$Temperature*1000 - 273.15
uadf_pre_og$species <- "uv"

adf <- rbind(cadf_pre_og,uadf_pre_og)
p3 <- ggplot( data = adf, aes(x=Temperature, y=norm_mrme, color=species)) +
  geom_point(aes(color=species), size = 1)+
  geom_line(aes(x = Temperature, y=fit))+
  geom_ribbon(data=adf, aes(x=Temperature, ymin=lwr, ymax=upr, fill=species), alpha=0.5, inherit.aes=F, show.legend = FALSE)+
  geom_vline(xintercept = 62.15, color = "#DA6C42FF", alpha = 0.5)+
  geom_vline(xintercept = 53.05, color = "#225BB2FF", alpha = 0.5)+
  scale_color_manual(labels = c("S. cerevisiae", "S. uvarum"),values=c('#DA6C42FF','#2166ac'))+
  ylab("Normalized secondary structure")+
  xlab(paste("Temperature ","\u00b0","C",sep=""))+
  ggtitle("AHA1") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title = element_blank())+
  theme(legend.position = "bottom")
p3


ggsave("Fig4.pdf",p1, width = 6, height = 4, units = "in")



##Wavelength spectra # supplementary
wv_acer <- read.csv("cAHA1_rep1_wv_edited.csv", header = T)
colnames(wv_acer) <- c("wavelength","4C","20C","30C","40C","50C","60C","75C")

wv_acer <- wv_acer %>%
  pivot_longer(!wavelength, names_to = "Temp", values_to = "Theta")
wv_acer$wavelength <- as.numeric(wv_acer$wavelength)
wv_acer$Theta <- as.numeric(wv_acer$Theta)

wv_acer$Temp <- factor(wv_acer$Temp, levels = c("4C",  "20C", "30C", "40C", "50C", "60C", "75C"))

#Calculate mean residue molar ellipticity
wv_acer$molar_ellip <- 100*wv_acer$Theta/(0.1*cA_M*1000)
wv_acer$mean_residue_molar_ellip <- wv_acer$molar_ellip/350


c2 <- ggplot( data = wv_acer, aes(x=wavelength, y=mean_residue_molar_ellip, group = Temp, color = Temp)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(n = 7, name = "Reds"),labels = c(paste("4","\u00b0","C", sep = ""),paste("20","\u00b0","C", sep = ""),paste("30","\u00b0","C", sep = ""),paste("40","\u00b0","C", sep = ""),paste("50","\u00b0","C", sep = ""),paste("60","\u00b0","C", sep = ""),paste("75","\u00b0","C", sep = "")))+
  ylab("Mean residue molar ellipticity")+
  xlab("Wavelength (nm)")+
  ggtitle("AHA1")+
  ylim(-9500,15)+
  theme(plot.title = element_text(hjust = 0.5))

c2


wv_auv <- read.csv("uAHA1_wv_melt_edited.csv", header = T)
colnames(wv_auv) <- c("wavelength","4C","20C","30C","40C","50C","60C","75C")

wv_auv <- wv_auv %>%
  pivot_longer(!wavelength, names_to = "Temp", values_to = "Theta")
wv_auv$wavelength <- as.numeric(wv_auv$wavelength)
wv_auv$Theta <- as.numeric(wv_auv$Theta)

wv_auv$Temp <- factor(wv_auv$Temp, levels = c("4C",  "20C", "30C", "40C", "50C", "60C", "75C"))

wv_auv <- wv_auv %>%
  filter(!wavelength < 210)

#Calculate mean residue molar ellipticity
wv_auv$molar_ellip <- 100*wv_auv$Theta/(0.1*uA_M*1000)
wv_auv$mean_residue_molar_ellip <- wv_auv$molar_ellip/347

u2 <- ggplot( data = wv_auv, aes(x=wavelength, y=mean_residue_molar_ellip, group = Temp, color = Temp)) +
  geom_line() +
  scale_color_manual(values = brewer.pal(n = 7, name = "Blues"),labels = c(paste("4","\u00b0","C", sep = ""),paste("20","\u00b0","C", sep = ""),paste("30","\u00b0","C", sep = ""),paste("40","\u00b0","C", sep = ""),paste("50","\u00b0","C", sep = ""),paste("60","\u00b0","C", sep = ""),paste("75","\u00b0","C", sep = "")))+
  ylab("Mean residue molar ellipticity")+
  xlab("Wavelength (nm)")+
  ggtitle("AHA1")+
  ylim(-9500,15)+
  theme(plot.title = element_text(hjust = 0.5))
u2

q2 <- ggarrange(c1,u1,c2,u2, ncol = 2, nrow = 2, widths=(c(1,1)), labels = c("A", "B", "C","D"))

q3 <- ggarrange(q2,p3, nrow = 2,labels = c("","E"), heights = (c(2,1)), widths = c(1,0.8))

ggsave("FigS6.pdf",q3, width = 8, height = 9, units = "in")

