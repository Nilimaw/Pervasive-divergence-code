##TPP analysis for Replicate 1, Replicate 2 and Hybrid. 
##To plot fitted curves from the three experiments.
##Normalization done using Tms for all genes.
##Hard code to analyze plots  to account for changes in curve fitting.
##Differences in R results in differing curve fitting

library (mstherm)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(reshape2)
sessionInfo()

###Analysis of TPP1
control <- system.file("extdata", "demo_project/nw_control_tpp1.tsv", package = "mstherm")
expt <- MSThermExperiment(control)

#Normalize to standard
par(mfrow=c(1,2))
expt1 <- normalize_to_std(expt, "CON__P02769")

res1 <- model_experiment(expt1,bootstrap = T, smooth = T, min_rep_psm = 3, np = 1)

result<- as.data.frame(res1)
result$cer1.tm.CI_l <- NA
result$cer1.tm.CI_u <- NA
result$uv1.tm.CI_l <- NA
result$uv1.tm.CI_u <- NA
result$gene <- rownames(result)


#Extract confidence intervals from dataset
for (g in result$gene){
  x<- res1[[g]]
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


write.csv(result,"C:/Users/Nilima/Desktop/mstherm_modified/mstherm_results_tpp1.csv", row.names = TRUE)

tpp1 <- read.csv("C:/Users/Nilima/Desktop/mstherm_modified/mstherm_results_tpp1.csv",header = T) 

###Analyze tpp2 data
control <- system.file("extdata", "demo_project/nw_control_tpp2.tsv", package = "mstherm")
expt <- MSThermExperiment(control)

#Normalize to standard
expt <- normalize_to_std(expt, "CON__P02769")
expt2 <- expt

##Print normed temps
expt2$samples$uv$replicates$uv2$meta$temp
expt2$samples$cer$replicates$cer2$meta$temp

##Use corrected temps
old <- c(30.0, 33.6, 37.7, 42.4, 45.7, 50.0, 54.8, 59.3, 64.8, 69.9)
new <- c(39.90096, 41.52294, 43.37019, 45.48778, 46.97459, 48.91196, 51.07460, 53.10207, 55.58010, 57.87790)
l <- lm(new ~ old)

s <- "cer2"
r <- "cer2"
expt2$samples[[s]]$replicates[[r]]$meta$temp <- expt2$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]


s <- "uv2"
r <- "uv2"
expt2$samples[[s]]$replicates[[r]]$meta$temp <- expt2$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]

expt2$samples$uv$replicates$uv2$meta$temp
expt2$samples$cer$replicates$cer2$meta$temp


#inter-replicate normalization - fit melting curves
data_norm2 <- model_experiment(expt2,bootstrap = T, smooth = T, min_rep_psm = 3, np = 1)
                                                                                                                       
norm_res <- as.data.frame(data_norm2)
norm_res$gene <- rownames(norm_res)

norm_res$cer2.tm.CI_l <- NA
norm_res$cer2.tm.CI_u <- NA
norm_res$uv2.tm.CI_l <- NA
norm_res$uv2.tm.CI_u <- NA


#Extract CIs for tpp2
for (g in norm_res$gene){
  x<- data_norm2[[g]]
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


tpp2n <- norm_res

write.csv(norm_res, "C:/Users/Nilima/Desktop/mstherm_modified/mstherm_tpp2_norm.csv", row.names = TRUE)


###Now normalizing hybrid subset1 <- uniquely mapping peptides
control <- system.file("extdata", "demo_project/nw_control_hy_subset1.tsv", package = "mstherm")
expt <- MSThermExperiment(control)

par(mfrow=c(1,2))
expt <- normalize_to_std(expt, "CON__P02769")

expt3 <- expt


##Print normed temps
expt3$samples$uv$replicates$uvh$meta$temp
expt3$samples$cer$replicates$cerh$meta$temp


#Use corrected temperatures
old <- c(30.0, 31.3, 33.6, 37.7, 42.4, 45.7, 48.3, 50.0, 52.0, 53.4, 55.9, 59.5, 64.0, 67.9, 70.5, 72.0)
new <- c(37.82327, 38.68449, 40.20817, 42.92431, 46.03794, 48.22410, 49.94653, 51.07273, 52.39768, 53.32514, 54.98132,57.36622, 60.34735, 62.93100, 64.65343, 65.64713)

l <- lm (new ~old)

s <- "cerh"
r <- "cerh"
expt3$samples[[s]]$replicates[[r]]$meta$temp <- expt3$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]


s <- "uvh"
r <- "uvh"
expt3$samples[[s]]$replicates[[r]]$meta$temp <- expt3$samples[[s]]$replicates[[r]]$meta$temp * l$coefficients[2] + l$coefficients[1]

expt3$samples$uv$replicates$uvh$meta$temp
expt3$samples$cer$replicates$cerh$meta$temp


#inter-replicate normalization - hybrid to tpp1
data_normh <- model_experiment(expt3,bootstrap = T, smooth = T, min_rep_psm = 3, np = 1)
norm_res <- as.data.frame(data_normh)
norm_res$gene <- rownames(norm_res)

norm_res$cerh.tm.CI_l <- NA
norm_res$cerh.tm.CI_u <- NA
norm_res$uvh.tm.CI_l <- NA
norm_res$uvh.tm.CI_u <- NA

#Extract CIs from hybrid data
for (g in norm_res$gene){
  x<- data_normh[[g]]
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

tpphn <- norm_res

write.csv(tpphn, "C:/Users/Nilima/Desktop/mstherm_modified/mstherm_tpph_norm.csv", row.names = TRUE)

setwd("C:/Users/Nilima/Desktop/mstherm_modified")


##Intersect the three datasets
a <- tpp1 %>%
  filter(uv1.r2 > 0.8 & cer1.r2 > 0.8)
b <- tpp2n %>%
  filter(cer2.r2 > 0.8 & uv2.r2 > 0.8 )
c <- tpphn %>%
  filter(cerh.r2 > 0.8 & uvh.r2 > 0.8)


tpp <- full_join(a, b , by = "gene")
tpp <- left_join(tpp, c, by = "gene")

tpp <- tpp %>%
  filter(cer1.r2 > 0.8 & uv1.r2 > 0.8 & cer2.r2 > 0.8 & uv2.r2 > 0.8) %>%
  mutate(tm1_diff = cer1.tm - uv1.tm, tm2_diff = cer2.tm - uv2.tm, tmh_diff = cerh.tm - uvh.tm)

tpp$q1[tpp$tm1_diff > 0] <- "+"
tpp$q1[tpp$tm1_diff < 0] <- "-"

tpp$q2[tpp$tm2_diff > 0] <- "+"
tpp$q2[tpp$tm2_diff < 0] <- "-"

tpp$qh[tpp$tmh_diff > 0] <- "+"
tpp$qh[tpp$tmh_diff < 0] <- "-"

tpp <- tpp %>%
  mutate(q12 = paste(q1,q2), q12h = paste(q12,qh))

table(tpp$q12)
table(tpp$q12h)

hy <- tpp %>%
  filter(!is.na(tmh_diff))

table(hy$q12h)


###Read list of proteins - main dataset
fil_data <- read.csv("proteomics_data.csv")
table(fil_data$Q1Q2)
hy_og <- fil_data %>%
  filter(!is.na(tmh_diff))

setdiff(fil_data$gene, tpp$gene)
#"YFR032C-A" "YGR034W"   "YNL004W"   "YOL040C"   "YDR381W"

missing <- c("YFR032C-A", "YGR034W",   "YNL004W",   "YOL040C",   "YDR381W")

setdiff(tpp$gene,fil_data$gene)
#"YLR028C" "YLR303W" "YLL018C" "YMR120C" "YBR026C" "YER156C"

gen_list <- c("YFR032C-A", "YGR034W", "YNL004W", "YOL040C", "YDR381W", "YLR028C", "YLR303W", "YLL018C", "YMR120C", "YBR026C", "YER156C")

setdiff(hy$gene, hy_og$gene)
#"YLR028C" "YLR303W" "YLL018C" "YMR120C" "YBR026C" "YER156C"

setdiff(hy_og$gene, hy$gene)
#"YNL004W" "YIL018W" "YOL040C" "YDR381W"

pdf("differences.pdf", 15, 5, pointsize=10)

for (gen in gen_list){
  print(gen)
  par(mfrow=c(1,3))
  plot(res1[[gen]])
  plot(data_norm2[[gen]])
  if (gen %in% tpphn$gene){
    plot(data_normh[[gen]])
  }
}
dev.off()

##Functions with modifications from mstherm
gen_profile <- function( x, method='sum', method.denom='first' ) {
  
  summarized <- switch( method,
                        sum    = apply(x, 2, sum),
                        median = apply(x, 2, median),
                        # for these, the inner apply() will transpose, so the outer is applied
                        # to rows rather than columns
                        ratio.median = 
                          apply( apply(x,1,abs_to_ratio,method=method.denom),1,median ),
                        ratio.mean = 
                          apply( apply(x,1,abs_to_ratio,method=method.denom),1,mean ),
                        stop("Invalid method type", call.=T)
  )
  return( abs_to_ratio(summarized, method=method.denom) )
  
}

abs_to_ratio <- function(x, method='first') {
  
  denom <- switch( method,
                   first  = x[1],
                   max    = max(x),
                   top3   = mean( x[ order(x,decreasing=T)[1:3] ] ),
                   near   = median( x[ x > x[1]*0.8 ] ),
                   #compat = {
                   #m <- mean( x[ order(x,decreasing=T)[1:3] ] )
                   #b <- mean( x[ x > m*0.8 & x < m*1.2 ] )
                   #ifelse(is.na(b),m,b)
                   #},
                   stop("Invalid method type",call.=T)
  )
  
  return( x/denom )
  
}

sigmoid <- function(p,k,m,x) {
  
  (1-p)/(1+exp(-k*(1/x-1/m)))+p
  
}

sigmoid.d1 <- function(p,k,m,x) {
  
  -((1 - p) * (exp(-k * (1/x - 1/m)) * (k * (1/x^2)))/(1 + exp(-k * (1/x - 1/m)))^2)
  
}

# look for probable missing values
is_consistent <- function(v,cutoff=0.3) {
  
  len <- length(v)
  
  if (len < 2) {
    return(1)
  }
  
  for (i in 1:(len-1)) {
    if (i == 1 & (v[i]/v[i+1]) < cutoff) {
      return(0)
    }
    if (i > 1) {
      if ((v[i]/v[i-1]) < cutoff & (v[i]/v[i+1]) < cutoff) {
        return(0)
      }
    }
  }
  
  return(1)
  
}


# Perform the actual model fitting
try_fit <- function(ratios,temps,trim,smooth) {
  
  x <- temps
  y <- ratios
  
  if (!missing(trim) & trim) {
    x <- temps[which.max(temps):length(temps)]
    y <- ratios[which.max(temps):length(temps)]
  }
  
  if (smooth) {
    f <- loess(y ~ x, span=0.65)
    y <- f$fitted
  }
  
  fit <- list()
  
  st.coarse <- expand.grid(p=c(0,0.3),k=seq(0,4000,by=1000),m=seq(30,60,by=15))
  st.fine   <- expand.grid(p=c(0,0.3),k=seq(0,8000,by=200),m=seq(30,80,by=10))
  for (st in list(st.coarse,st.fine)) {
    tryCatch( {
      mod <- nls2::nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=st,algorithm="brute-force",control=nls.control(warnOnly=T,maxiter=5000))
      fit <-
        nls2::nls2(y~sigmoid(p,k,m,x),data=list(x=x,y=y),start=mod,control=nls.control(warnOnly=T),algorithm="port",lower=c(0,1,10),upper=c(0.4,100000,100))
      obj <- list()
      obj$plat  <- as.numeric(coefficients(fit)[1])
      obj$k     <- as.numeric(coefficients(fit)[2])
      obj$tm    <- as.numeric(coefficients(fit)[3])
      obj$slope <- as.numeric(sigmoid.d1(obj$plat,obj$k,obj$tm,obj$tm))
      y.fit <- sigmoid(obj$plat,obj$k,obj$tm,temps)
      obj$y.fit <- y.fit
      obj$resid <- ratios - y.fit
      obj$r2 <- 1-(sum(obj$resid^2)/(length(ratios)*var(ratios)))
      obj$rmsd <- sqrt( sum(obj$resid^2)/length(ratios) )
      return(obj)
    },error = function(e) {})
  }
  return(NULL)
  
}

model_protein <- function( expt, protein,min_rep_psm  = 0,
                           min_smp_psm  = 0,
                           min_tot_psm  = 0,
                           max_inf      = 1,
                           min_score,
                           max_score,
                           smooth       = 0,
                           method       = 'sum',
                           method.denom = 'near',
                           trim         = 0,
                           bootstrap    = 0,
                           min_bs_psms  = 8,
                           annot_sep    = '|',
                           max_slope    = 0,
                           min_r2       = 0,
                           min_reps     = 0,
                           only_modeled = 0,
                           check_missing = 0,
                           missing_cutoff = 0.3
) {
  
  self <- structure(
    list(
      name         = protein,
      series       = list(),
      sample_names = c(),
      parameters   = as.list(environment())
    ),
    class = "MSThermResult"
  )
  self$parameters[['expt']] <- NULL
  
  #self$annotation <- gen_description(expt, protein, sep=annot_sep)
  
  psm_tot   <- 0 # track total PSMs for protein
  n_samples <- length(expt$samples)
  
  if (min_r2 > 0 || max_slope < 0 ) {
    only_modeled = 1
  }
  
  for (i_sample in 1:n_samples) {
    
    psm_smp <- 0 # track total PSMs for sample
    worst_slope <- -1
    worst_r2    <- 1
    
    sample <- expt$samples[[i_sample]]
    n_replicates <- length(sample$replicates)
    self$sample_names[i_sample] <- sample$name
    
    n_reps <- 0
    
    # Process each replicate
    for (i_replicate in 1:n_replicates) {
      
      replicate <- sample$replicates[[i_replicate]]
      
      # Set score range if not yet defined
      if ("score" %in% colnames(replicate$data)) {
        if (missing(max_score)) {
          max_score = max(replicate$data$score)
        }
        if (missing(min_score)) {
          min_score = min(replicate$data$score)
        }
      }
      
      temps <- replicate$meta$temp
      
      # Track global min and max temperature points for the protein
      if (is.null(self$tmin) || self$tmin > min(temps)) {
        self$tmin <- min(temps)
      }
      if (is.null(self$tmax) || self$tmax < max(temps)) {
        self$tmax <- max(temps)
      }
      
      # Pull out matching data points passing coisolation and score
      # thresholds
      sub <- replicate$data[which(replicate$data$protein == protein
                                  & replicate$data$coelute_inf <= max_inf),];
      if ("score" %in% colnames(sub)) {
        sub <- sub[which( sub$score >= min_score & sub$score <= max_score ),]
      }
      
      if (nrow(sub) < 1) {
        next
      }
      
      # Reorder channels based on metadata
      quant_columns <- match(replicate$meta$channel,colnames(sub))
      quant <- sub[,quant_columns]
      
      # Filter out rows with NA or with all zero values
      ok <- apply(quant,1,function(v) {
        all(!is.na(v)) &
          any(v>0)     &
          ( (! check_missing) || is_consistent(v, missing_cutoff) )
      })
      sub   <- sub[ok,]
      quant <- quant[ok,]
      
      n_psms  <- nrow(sub)
      
      # Update PSM totals and check cutoffs
      psm_tot <- psm_tot + n_psms
      psm_smp <- psm_smp + n_psms
      if (n_psms < min_rep_psm) {
        return(NULL)
      }
      
      # Obviously don't try to model if no rows pass filtering
      if (n_psms < 1) {
        next
      }
      
      profile <- gen_profile(quant,method,method.denom=method.denom)
      
      fit <- try_fit(profile, temps, trim=trim, smooth=smooth)
      fit$is.fitted <- !is.null(fit)
      
      # calculate weighted co-inf
      sums       <- apply(quant, 1, sum)
      fit$inf    <- sum(sub$coelute_inf * sums) / sum(sums)
      
      # keep track of other data for later use
      fit$psm    <- n_psms
      fit$name   <- replicate$name
      fit$sample <- sample$name
      fit$x      <- temps
      fit$y      <- profile
      
      if (fit$is.fitted) {
        
        #bootstrap if asked
        bs <- c()
        iterations <- 20
        bs.ratios <- matrix(nrow=iterations,ncol=length(profile))
        fit.count <- 0
        if (bootstrap & nrow(quant) >= min_bs_psms) { 
          for (i in 1:iterations) {
            
            quant.bs <- quant[sample(nrow(quant),nrow(quant),replace=T),]
            profile.bs <- gen_profile(quant.bs,method,method.denom=method.denom)
            bs.ratios[i,] <- profile.bs
            fit.bs <- try_fit(profile.bs,temps,trim,smooth)
            is.fitted <- !is.null(fit.bs)
            if (is.fitted) {
              bs[i] <- fit.bs$tm
              fit.count <- fit.count + 1
            }
          }
          fit$bs.lowers <- apply(bs.ratios,2,function(x) quantile(x,0.025))
          fit$bs.uppers <- apply(bs.ratios,2,function(x) quantile(x,0.975))
          if (fit.count > iterations * 0.8) {
            fit$tm_CI <- quantile(bs,c(0.025,0.975),na.rm=T)
          }
        }
        if (fit$r2 < min_r2) {
          next;
        }
        if (fit$slope > max_slope) {
          next;
        }
        n_reps = n_reps + 1
        
      }
      
      # If unfitted, these variables still need to be set but should be NA
      else {
        if (only_modeled) {
          next;
        }
        fit$tm    <- NA
        fit$slope <- NA
        fit$k     <- NA
        fit$plat  <- NA
        fit$r2    <- NA
      }
      
      self$series[[replicate$name]] <- fit
    }
    
    # Filter on total PSMs for sample
    if (psm_smp < min_smp_psm) {
      return( NULL )
    }
    # Filter on minimum modeled replicates for sample
    if (n_reps < min_reps) {
      return( NULL )
    }
    
  }
  
  # Do final filtering on total PSMs for protein
  if (psm_tot < min_tot_psm) {
    return( NULL )
  }
  # Filter on worst slope
  if (worst_slope > max_slope) {
    return( NULL )
  }
  # Filter on worst R2
  if (worst_r2 < min_r2) {
    return( NULL )
  }
  
  return( self )
}



##Force fit missing proteins
model_YFR032CA   <- model_protein(expt1, "YFR032C-A", smooth=TRUE, bootstrap=TRUE)
model_YGR034W   <- model_protein(expt1, "YGR034W", smooth=TRUE, bootstrap = TRUE)
model_YNL004W   <- model_protein(expt1, "YNL004W", smooth=TRUE, bootstrap=TRUE)
model_YOL040C   <- model_protein(expt1, "YOL040C", smooth=TRUE, bootstrap=TRUE)
model_YDR381W   <- model_protein(expt1, "YDR381W", smooth=TRUE, bootstrap=TRUE)

YFR032CA <- unlist(model_YFR032CA)
YGR034W <- unlist(model_YGR034W)
YNL004W <- unlist(model_YNL004W)
YOL040C <- unlist(model_YOL040C)
YDR381W <- unlist(model_YDR381W)

##Missing in hybrid
modelh_YIL018W <- model_protein(expt3, "YIL018W", smooth=TRUE, bootstrap=TRUE)


pdf("filtered_plots_all.pdf", 15, 5, pointsize=10)
for (xgen in fil_data$gene){
  print(xgen)
  if(xgen %in% missing){
    next
  }
  if(xgen == "YIL018W"){
    next
  }
  if (!is.na(fil_data$tmh_diff[fil_data$gene == xgen])){
    par(mfrow=c(1,3))
    plot(res1[[xgen]])
    plot(data_norm2[[xgen]])
    plot(data_normh[[xgen]])
  }
  else{
    par(mfrow = c(1,2))
    plot(res1[[xgen]])
    plot(data_norm2[[xgen]])
  }
}


#YFR032C-A
par(mfrow=c(1,2))
xgen <- "YFR032C-A"
plot(model_YFR032CA)
plot(data_norm2[[xgen]])

#YGR034W
par(mfrow=c(1,2))
xgen <- "YGR034W"
plot(model_YGR034W)
plot(data_norm2[[xgen]])

#YNL004W  
par(mfrow=c(1,2))
xgen <- "YNL004W"
plot(model_YNL004W)
plot(data_norm2[[xgen]])


#YOL040C
par(mfrow=c(1,3))
xgen <- "YOL040C"
plot(model_YOL040C)
plot(data_norm2[[xgen]])
plot(data_normh[[xgen]])

#YDR381W
par(mfrow=c(1,3))
xgen <- "YDR381W"
plot(model_YDR381W)
plot(data_norm2[[xgen]])
plot(data_normh[[xgen]])


##YIL018W
par(mfrow=c(1,3))
xgen <- "YIL018W"
plot(res1[[xgen]])
plot(data_norm2[[xgen]])
plot(modelh_YIL018W)


dev.off()


###Plot Tm fits for proteins whose enzyme activity was assayed

#MDH1: YKL085W
#MDH3: YDL078C
#GLR1: YPL091W
#HXK1: YFR053C
#HXK2: YGL253W
#GLK1: YCL040W
#EMI2: YDR516C

enz <- c("YKL085W", "YDL078C", "YPL091W", "YFR053C", "YGL253W","YCL040W","YDR516C")
proteins <- c("MDH1","MDH3","GLR1","HXK1","HXK2","GLK1","EMI2")

##Plot FigS1
pdf("FigS1.pdf",10,25,pointsize = 10)
par(mfrow=c(7,2))
for (xgen in enz){
  print(xgen)
  plot(res1[[xgen]])
  plot(data_norm2[[xgen]])
}

dev.off()

