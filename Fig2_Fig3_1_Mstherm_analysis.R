#TPP analysis for Replicate 1, Replicate 2 and Hybrid. 
##Normalization done using Tms for all genes

library (mstherm)
library(plyr)
library(ggplot2)

###Analysis of TPP1
control <- system.file("extdata", "demo_project/nw_control_tpp1.tsv", package = "mstherm")
expt <- MSThermExperiment(control)

#Normalize to standard
par(mfrow=c(1,2))
expt <- normalize_to_std(expt, "CON__P02769")
res <- model_experiment(expt,bootstrap = T, smooth = T, min_rep_psm = 3, np = 1)

result<- as.data.frame(res)
result$cer1.tm.CI_l <- NA
result$cer1.tm.CI_u <- NA
result$uv1.tm.CI_l <- NA
result$uv1.tm.CI_u <- NA
result$gene <- rownames(result)


#Extract confidence intervals from dataset
for (g in result$gene)
{
  x<- res[[g]]
  if (!is.null(x$series$cer1$tm_CI))
  {
    result[g,"cer1.tm.CI_l"] <- x$series$cer1$tm_CI[1]
    result[g,"cer1.tm.CI_u"] <- x$series$cer1$tm_CI[2]
  }
  if (!is.null(x$series$uv1$tm_CI))
  {
    result[g,"uv1.tm.CI_l"] <- x$series$uv1$tm_CI[1]
    result[g,"uv1.tm.CI_u"] <- x$series$uv1$tm_CI[2]
  }
  
}

write.csv(result,"C:/Users/nwalu/Documents/Fay_lab_work_2020/Aim2_Mass spec/2021_August_Mass_spec_analysis with cleaned ortholog set/21Dec_TPP12H/all/mstherm_results_tpp1.csv", row.names = TRUE)

##Set tpp1 as result df
tpp1 <- read.csv("C:/Users/nwalu/Documents/Fay_lab_work_2020/Aim2_Mass spec/2021_August_Mass_spec_analysis with cleaned ortholog set/21Dec_TPP12H/all/mstherm_results_tpp1.csv",header = T) 


###Analyze tpp2 data
control <- system.file("extdata", "demo_project/nw_control_tpp2.tsv", package = "mstherm")
expt <- MSThermExperiment(control)

#Normalize to standard
par(mar=c(1,1,1,1))
par(mfrow=c(2,3))
expt1 <- normalize_to_std(expt, "CON__P02769", plot = FALSE)
res <- model_experiment(expt1, smooth=TRUE, bootstrap=FALSE, min_rep_psm=3, np=2)

expt2 <- expt1

tpp2 <- as.data.frame(res)
tpp2$X <- rownames(tpp2)

##Cer Tms
tpp1_c <- tpp1[c("X", "cer1.tm","cer1.r2")]
tpp2_c <- tpp2[c("X", "cer2.tm","cer2.r2")]

tpp12_c <- merge(tpp1_c,tpp2_c,by = "X")

##Uv tms
tpp1_u <- tpp1[c("X", "uv1.tm","uv1.r2")]
tpp2_u <- tpp2[c("X", "uv2.tm","uv2.r2")]

tpp12_u <- merge(tpp1_u,tpp2_u,by = "X")


##Baseline Tms from combining cer1 and uv1
tm.b1 <- tpp12_c$cer1.tm
tm.b2 <- tpp12_u$uv1.tm

tm.b1 <- as.data.frame(unlist(tm.b1))
tm.b2 <- as.data.frame(unlist(tm.b2))
colnames(tm.b1) <- "tm"
colnames(tm.b2) <- "tm"
tm.b <- rbind(tm.b1,tm.b2)
tm.b <- unlist(tm.b)


##Fit tms for tpp2 (cer2 and uv2) to tpp1
tm.u1 <- tpp12_c$cer2.tm
tm.u2 <- tpp12_u$uv2.tm

tm.u1 <- as.data.frame(unlist(tm.u1))
tm.u2 <- as.data.frame(unlist(tm.u2))
colnames(tm.u1) <- "tm"
colnames(tm.u2) <- "tm"
tm.u <- rbind(tm.u1,tm.u2)
tm.u <- unlist(tm.u)


##Fit linear model
l <- lm(tm.b ~ tm.u)

#Normalize temperatures
s <- "cer2"
r <- "cer2"
expt2$samples[[s]]$replicates[[r]]$meta$temp <- expt1$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]


s <- "uv2"
r <- "uv2"
expt2$samples[[s]]$replicates[[r]]$meta$temp <- expt1$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]

##Print normed temps
expt2$samples$uv$replicates$uv2$meta$temp
expt2$samples$cer$replicates$cer2$meta$temp
#39.90096 41.52294 43.37019 45.48778 46.97459 48.91196 51.07460 53.10207 55.58010 57.87790

#inter-replicate normalization - fit melting curves
data_norm <- model_experiment(expt2,bootstrap = T, smooth = T, min_rep_psm = 3, np = 1)
norm_res <- as.data.frame(data_norm)
norm_res$gene <- rownames(norm_res)

norm_res$cer2.tm.CI_l <- NA
norm_res$cer2.tm.CI_u <- NA
norm_res$uv2.tm.CI_l <- NA
norm_res$uv2.tm.CI_u <- NA


#Extract CIs for tpp2
for (g in norm_res$gene)
{
  x<- data_norm[[g]]
  if (!is.null(x$series$uv2$tm_CI))
  {
    norm_res[g,"uv2.tm.CI_l"] <- x$series$uv2$tm_CI[1]
    norm_res[g,"uv2.tm.CI_u"] <- x$series$uv2$tm_CI[2]
  }
  if (!is.null(x$series$cer2$tm_CI))
  {
    norm_res[g,"cer2.tm.CI_l"] <- x$series$cer2$tm_CI[1]
    norm_res[g,"cer2.tm.CI_u"] <- x$series$cer2$tm_CI[2]
  }
  
}


write.csv(norm_res, "C:/Users/nwalu/Documents/Fay_lab_work_2020/Aim2_Mass spec/2021_August_Mass_spec_analysis with cleaned ortholog set/21Dec_TPP12H/all/mstherm_tpp2_norm.csv", row.names = TRUE)


###Now normalizing hybrid subset1 <- uniquely mapping peptides
control <- system.file("extdata", "demo_project/nw_control_hy_subset1.tsv", package = "mstherm")
expt <- MSThermExperiment(control)

par(mfrow=c(1,2))
expt1 <- normalize_to_std(expt, "CON__P02769")
res <- model_experiment(expt1,bootstrap = F, smooth = T, min_rep_psm = 3, np = 1)

expt2 <- expt1

tpph <- as.data.frame(res)
tpph$X <- rownames(tpph)

##Cer Tms
tpp1_c <- tpp1[c("X", "cer1.tm","cer1.r2")]
tpph_c <- tpph[c("X", "cerh.tm","cerh.r2")]

tpp1h_c <- merge(tpp1_c,tpph_c,by = "X")

##Uv tms
tpp1_u <- tpp1[c("X", "uv1.tm","uv1.r2")]
tpph_u <- tpph[c("X", "uvh.tm","uvh.r2")]

tpp1h_u <- merge(tpp1_u,tpph_u,by = "X")

##Baseline Tms from combining cer1 and uv1
tm.b1 <- tpp1h_c$cer1.tm
tm.b2 <- tpp1h_u$uv1.tm

tm.b1 <- as.data.frame(unlist(tm.b1))
tm.b2 <- as.data.frame(unlist(tm.b2))
colnames(tm.b1) <- "tm"
colnames(tm.b2) <- "tm"
tm.b <- rbind(tm.b1,tm.b2)
temp_1 <- tm.b
tm.b <- unlist(tm.b)

##Fit tms for tpp2 cer2 and uv2
tm.u1 <- tpp1h_c$cerh.tm
tm.u2 <- tpp1h_u$uvh.tm

tm.u1 <- as.data.frame(unlist(tm.u1))
tm.u2 <- as.data.frame(unlist(tm.u2))
colnames(tm.u1) <- "tm"
colnames(tm.u2) <- "tm"
tm.u <- rbind(tm.u1,tm.u2)
temp_h <- tm.u
tm.u <- unlist(tm.u)

colnames(temp_1)[1] <- "Tm_rep1"
temp_1$species <- "cer"
temp_1$species[691:1382] <- "uv"

colnames(temp_h)[1] <- "Tm_reph"
temp_h$species <- "cer"
temp_h$species[691:1382] <- "uv"

df <- cbind(temp_1,temp_h)
df <- df[,1:3]

p <- ggplot(df, aes(x=Tm_reph, y=Tm_rep1 )) +
  geom_point(aes(color = species))+
  geom_smooth(method = 'lm')+
  ggtitle("Pre_norm_Tms")+
  labs(x = "TPP1 tms", y = "TPPH - Tms" )
p
ggsave("Fay_lab_work_2020/Aim2_Mass spec/2021_August_Mass_spec_analysis with cleaned ortholog set/21Dec_TPP12H/all/Pre_norm_tms.png",p)


##Fit linear model
l <- lm(tm.b ~ tm.u)

s <- "cerh"
r <- "cerh"
expt2$samples[[s]]$replicates[[r]]$meta$temp <- expt1$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]


s <- "uvh"
r <- "uvh"
expt2$samples[[s]]$replicates[[r]]$meta$temp <- expt1$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]

##Print normed temps
expt2$samples$uv$replicates$uvh$meta$temp
expt2$samples$cer$replicates$cerh$meta$temp

##37.82327 38.68449 40.20817 42.92431 46.03794 48.22410 49.94653 51.07273 52.39768 53.32514 54.98132
## 57.36622 60.34735 62.93100 64.65343 65.64713

#inter-replicate normalization - hybrid to tpp1
data_norm <- model_experiment(expt2,bootstrap = T, smooth = T, min_rep_psm = 3, np = 1)
norm_res <- as.data.frame(data_norm)
norm_res$gene <- rownames(norm_res)

norm_res$cerh.tm.CI_l <- NA
norm_res$cerh.tm.CI_u <- NA
norm_res$uvh.tm.CI_l <- NA
norm_res$uvh.tm.CI_u <- NA

#Extract CIs from hybrid data
for (g in norm_res$gene)
{
  x<- data_norm[[g]]
  if (!is.null(x$series$uvh$tm_CI))
  {
    norm_res[g,"uvh.tm.CI_l"] <- x$series$uvh$tm_CI[1]
    norm_res[g,"uvh.tm.CI_u"] <- x$series$uvh$tm_CI[2]
  }
  if (!is.null(x$series$cerh$tm_CI))
  {
    norm_res[g,"cerh.tm.CI_l"] <- x$series$cerh$tm_CI[1]
    norm_res[g,"cerh.tm.CI_u"] <- x$series$cerh$tm_CI[2]
  }
  
}

write.csv(norm_res, "C:/Users/nwalu/Documents/Fay_lab_work_2020/Aim2_Mass spec/2021_August_Mass_spec_analysis with cleaned ortholog set/21Dec_TPP12H/all/mstherm_tpph_norm.csv", row.names = TRUE)

