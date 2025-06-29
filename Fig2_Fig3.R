###TPP12H normed plots
###Figure 2 and Figure 3 new

setwd("C:/Users/Nilima/Desktop/mstherm_modified")

library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(showtext)
library(patchwork)
library(GGally)
library(dplyr)
library(gridExtra)
library(hrbrthemes)
library(ggpubr)
library(stringr)
library(aod)


##Read data with added annotations
loc1 <- "proteomics_data.csv"
## Same as final_dataset_TPP12H_normed.csv 

data <- read_csv(loc1)
colnames(data)
data <- data%>%
  mutate(avg_cer_tm = (cer1.tm + cer2.tm)/2) %>%
  mutate(avg_uv_tm = (uv1.tm + uv2.tm)/2) %>%
  mutate(tm12_diff_avg = (tm1_diff + tm2_diff)/2)

data$TPP1_int_TPP2_NOCI[data$TPP1_int_TPP2_NOCI == 1] <- "Non-overlapping confidence intervals"
data$TPP1_int_TPP2_NOCI[data$TPP1_int_TPP2_NOCI == 0] <- "NA"



theme_set(theme_minimal())
theme_update(
  # The size of the axes labels are different for x and y.
  axis.text.x = element_text(size = 6, margin = margin(t = 5)),
  axis.text.y = element_text(size = 6, margin = margin(r = 5)),
  axis.title.x = element_text(size = 8),
  axis.title.y = element_text(size = 8),
  plot.title = element_text(size = 6, face = "bold",margin = margin(t = 15)),
  legend.text=element_text(size=15))


#Figure2

##Plot TPP1 vs TPP2
xlab <- "Tm difference (S.cer - S. uv, replicate 1)"
ylab <- "Tm difference (S.cer - S. uv, replictae 2)"

p1 <- ggplot(data, aes(x= tm1_diff, y=tm2_diff,color = TPP1_int_TPP2_NOCI)) +
  geom_point(size = 1)+
  scale_color_manual(name = "Confidence intervals", values=c("black", "orange"), breaks=c("", "Non-overlapping confidence intervals"), labels=c("","Non-overlapping confidence intervals" ))+
  labs(x = xlab, y = ylab )+
  geom_abline(aes(intercept = 0,slope = 0))+
  geom_vline(aes(xintercept = 0))+
  theme(legend.position = "none")+
  theme_bw()

p1
ggsave("Fig2.pdf",p1, width = 8, height = 4, units = "in")


###Figure 3: Histogram of hybrid - parent
##A: Tm difference in Parents vs Hybrid
data <- data %>% 
  mutate(hyminuscer = cerh.tm - avg_cer_tm) %>%
  mutate(hyminusuv = uvh.tm - avg_uv_tm)

theme_set(theme_classic())
theme_update(
  # The size of the axes labels are different for x and y.
  axis.text.x = element_text(size = 10, margin = margin(t = 5)),
  axis.text.y = element_text(size = 10, margin = margin(r = 5)),
  axis.title.x = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  plot.title = element_text(size = 10, face = "bold",margin = margin(t = 15)),
  legend.text=element_text(size=10))

##Subet dataset with proteins in rep1, rep2 and hybrid
hy <- data %>%
  filter(!is.na(tmh_diff ))

hy$tmpdiff_tmhdiff <- hy$tm12_diff - hy$tmh_diff

##Subset average Tm difference in parents
pl <- hy %>%
  select(tm12_diff)
pl$sp <- "Parents"
colnames(pl)[1]<- "tmdiff"

pl1 <- hy %>%
  select(tmh_diff)
pl1$sp <- "Hybrid"
colnames(pl1)[1] <- "tmdiff"

lp <- rbind(pl,pl1)

z1 <- lp %>%
  ggplot( aes(x=tmdiff, fill = sp)) +
  geom_histogram( color="grey39", linewidth = 0.3, alpha=0.6, position = 'identity') +
  geom_vline(xintercept = 1.8, color = "#c57f2f")+
  geom_vline(xintercept = 0.68, color = "#472d73")+
  scale_fill_manual(values=c("#5e3c99", "#f1a340")) +
  labs(fill="")+
  xlab("Tm difference (S. cer - S.uv)")+
  ylab("Count")+
  theme(legend.key.spacing.y = unit(0.4, 'cm'))

z1

wilcox.test(pl$tmdiff, pl1$tmdiff ,paired = TRUE, alternative="two.sided")$p.value



# Plot Figure 3B: HyTm - Parent Tm distributions for S. cer and S. uv
data_v3 <- as.data.frame(data$hyminuscer)
data_v3$species <- "cer"
colnames(data_v3)[1] <- "hyminus"
data_v3 <- data_v3 %>% 
  filter(!is.na(hyminus))
  

data_v4 <- as.data.frame(data$hyminusuv)
data_v4$species <- "uv"
colnames(data_v4)[1] <- "hyminus"
data_v4 <- data_v4 %>% 
  filter(!is.na(hyminus))


pl_data <- rbind(data_v3,data_v4)

z <- pl_data %>%
  ggplot( aes(x=hyminus, fill=species)) +
  geom_histogram( color="grey39", linewidth = 0.3, alpha=0.6, position = 'identity') +
  scale_fill_manual(labels = c("S. cerevisiae", "S. uvarum"), values=c("#DA6C42FF", "#225BB2FF")) +
  geom_vline(xintercept = -0.088, color = "#DA6C42FF")+
  geom_vline(xintercept = 1.032, color = "#225BB2FF")+
  labs(fill="")+
  ylab("Count")+
  xlab("Tm in hybrid - Tm in parent")+
  theme(legend.key.spacing.y = unit(0.4, 'cm'))
  
z

f2 <- ggarrange(z1, z,ncol = 2, nrow = 1,labels = c("A", "B"))
ggsave("Fig2_v4.pdf",f2, width = 10, height = 3, units = "in", bg = "white")


hy <- data %>%
  filter(qsummary == "+++" | qsummary == "---")
binom.test(179, 204, p = 0.5, alternative = "two.sided",conf.level = 0.95)$p.value


temp <- data %>%
  select(hyminuscer, hyminusuv, gene) %>%
  filter(!is.na(hyminuscer)) %>%
  filter(!is.na(hyminusuv))
t.test(temp$hyminuscer, temp$hyminusuv, paired = TRUE, alternative = "two.sided")$p.value



table(data$ribo)
table(data$ribo_GO)

r <- data %>%
  filter(q12 == "++" | q12 == "--")
table(r$ribo_GO)

table(r$ribo_GO,r$q12)
dat <- table(r$ribo_GO, r$q12)
fisher.test(dat)$p.value



data <- data %>%
  mutate(avgtm= (cer1.tm + uv1.tm + cer2.tm + uv2.tm)/4) %>%
  mutate(avgtmdiff = (tm1_diff + tm2_diff)/2)
t1 <- data$avgtmdiff
t2 <- data$avgtm

cor.test(t2, t1, method = "spearman") 
