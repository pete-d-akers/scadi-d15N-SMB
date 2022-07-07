#========================East Antarctic d15NNO3-SMB Analysis Script==========================
# This script analyzes and plots the surface mass balance (SMB) relationship with nitrate isotopes (d15N) 
#   at sites across East Antarctica, and reconstructs past SMBs from d15N in ice cores.
# Written by Pete D Akers, May 2019-May 2022

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(cowplot)
library(tidyverse)

#### Additional Functions ####
#Geom Stepribbon https://github.com/adibender/pammtools/blob/master/R/ggplot-extensions.R
geom_stepribbon <- function(
  mapping     = NULL,
  data        = NULL,
  stat        = "identity",
  position    = "identity",
  na.rm       = FALSE,
  show.legend = NA,
  inherit.aes = TRUE, ...) {
  
  layer(
    data        = data,
    mapping     = mapping,
    stat        = stat,
    geom        = GeomStepribbon,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, ... )
  )
  
}

#' @rdname geom_stepribbon
#' @format NULL
#' @usage NULL
#' @export

GeomStepribbon <- ggproto(
  "GeomStepribbon", GeomRibbon,
  
  extra_params = c("na.rm"),
  
  draw_group = function(data, panel_scales, coord, na.rm = FALSE) {
    
    if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
    data   <- rbind(data, data)
    data   <- data[order(data$x), ]
    data$x <- c(data$x[2:nrow(data)], NA)
    data   <- data[complete.cases(data["x"]), ]
    GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)
    
  }
  
)

lighten <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col + (255 - col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

darken <- function(color, factor = 0.5) {
  if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
  col <- col2rgb(color)
  col <- col - (col)*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

#### END Additional Functions ####

#### Loading nitrate data ####

# **FIRST DATA STEPS** 
  #Download file from https://doi.pangaea.de/10.1594/PANGAEA.941480 and rename to "east_antarctic_d15NNO3.tab" and delete 186 header lines
  #Download file from https://doi.pangaea.de/10.1594/PANGAEA.941490 and rename to "ABN1314-103_NO3_iso.tab" and delete 26 header rows
  #Download file from https://doi.pangaea.de/10.1594/PANGAEA.941489 and rename to "ABN1314-103_density.tab" and delete 23 header rows
all.nitrate <- read.table("east_antarctic_d15NNO3.tab", header=TRUE, sep="\t")
colnames(all.nitrate) <- c("event", "site", "alt.site", "latitude", "longitude", "elev.rema", "smb.ground",
                           "smb.mar.orig", "smb.mar.adj", "smb.mix", "type", "campaign", "transect", "cite.iso",
                           "cite.smb","smb.source", "NO3", "d15N", "note")
rownames(all.nitrate) <- NULL #resetting row numbers
all.nitrate$transect <- factor(all.nitrate$transect, levels=c("EAIIST", "DDU-Vostok", "JARE", "Zhongshan-DA", "CHICTABA", "Other")) #Setting desired order of transects
all.nitrate.nophoto <- all.nitrate[which(all.nitrate$d15N.nophoto == 1), ] #samples whose d15N are not different than SKL values (< +20)
all.nitrate.photo <- all.nitrate[-which(all.nitrate$d15N.nophoto == 1), ]
all.nitrate$NO3.flux <- all.nitrate$smb.mar.adj/all.nitrate$NO3
abn <- read.table("ABN1314-103_NO3_iso.tab", header=TRUE, sep="\t") #ABN201314 nitrate ice core data
colnames(abn) <- c("id", "depth.top", "depth.btm", "depth.mean", "yr.top", "yr.btm", "yr.mean", "yr.ka.BP",
                   "yr.int", "d15N", "d15N.sd", "NO3", "smb.gpr.upstream")
abn$yr.top.BP <- 1950-abn$yr.top
abn.density.all <- read.table("ABN1314-103_density.tab", header=TRUE, sep="\t") #ABN SMB from ice density
colnames(abn.density.all) <- c("id", "depth.top", "depth.btm", "depth.mean", "yr.top", "yr.btm", "yr.ka.BP", "yr.int", "density", "smb.density")
abn.density.all$yr.top.BP <- 1950-abn.density.all$yr.top
abn.density <- abn.density.all[-seq(1,3),] #Removing top layers with uncertain density

#Setting color palette for plotting
transect.palette <- brewer.pal(7, "Dark2")
transect.palette <- transect.palette[c(4,3,1,2,6)]
transect.palette[4] <- brewer.pal(7,"Set1")[7]
transect.palette[6] <- "dodgerblue"

#### END Loading Data ####

#####=====Calculating regressions #####
##All sites d15N-SMB Relationship
#Best mix (ground + Adj MAR for sites without ground)
lm.smb.d15N.mix <- lm(log(d15N+1)~I(1/(smb.mix)), data=all.nitrate) #I() is AsIs to make function read correctly
lm.smb.predict.d15N.mix <- lm(I(1/smb.mix)~log(d15N+1), data=all.nitrate)
lm.mix.sum <- summary(lm.smb.d15N.mix)
lm.mix.predict.sum <- summary(lm.smb.predict.d15N.mix)
smb.range.mix <- c((exp(coef(lm.smb.d15N.mix)[1])-1), #d15NNO3o
                   (exp(coef(lm.smb.d15N.mix)[2]/500 + coef(lm.smb.d15N.mix)[1])-1), #d15NNO3 at 500 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.mix)[2]/200 + coef(lm.smb.d15N.mix)[1])-1), #d15NNO3 at 200 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.mix)[2]/100 + coef(lm.smb.d15N.mix)[1])-1), #d15NNO3 at 100 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.mix)[2]/50 + coef(lm.smb.d15N.mix)[1])-1), #d15NNO3 at 50 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.mix)[2]/20 + coef(lm.smb.d15N.mix)[1])-1), #d15NNO3 at 20 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.mix)[2]/10 + coef(lm.smb.d15N.mix)[1])-1)) #d15NNO3 at 10 kg m-2 a-1
lm.coef.mix <- c(round(c(lm.mix.sum$coef[2],lm.mix.sum$coef[4],lm.mix.sum$coef[1],lm.mix.sum$coef[3],
                         lm.mix.sum$adj.r.squared), 2),
                 round(c(lm.mix.predict.sum$coef[2],lm.mix.predict.sum$coef[4],lm.mix.predict.sum$coef[1],
                         lm.mix.predict.sum$coef[3], lm.mix.predict.sum$adj.r.squared), 3))

##Comparison for alternate SMB source data
#Adjusted MAR for all sites
lm.smb.d15N.adj <- lm(log(d15N+1)~I(1/smb.mar.adj), data=all.nitrate)
lm.smb.predict.d15N.adj <- lm(I(1/smb.mar.adj)~log(d15N+1), data=all.nitrate)
lm.adj.sum <- summary(lm.smb.d15N.adj)
lm.adj.predict.sum <- summary(lm.smb.predict.d15N.adj)
smb.range.adj <- c((exp(coef(lm.smb.d15N.adj)[1])-1), #d15NNO3o
                   (exp(coef(lm.smb.d15N.adj)[2]/500 + coef(lm.smb.d15N.adj)[1])-1), #d15NNO3 at 500 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.adj)[2]/200 + coef(lm.smb.d15N.adj)[1])-1), #d15NNO3 at 200 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.adj)[2]/100 + coef(lm.smb.d15N.adj)[1])-1), #d15NNO3 at 100 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.adj)[2]/50 + coef(lm.smb.d15N.adj)[1])-1), #d15NNO3 at 50 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.adj)[2]/20 + coef(lm.smb.d15N.adj)[1])-1), #d15NNO3 at 20 kg m-2 a-1
                   (exp(coef(lm.smb.d15N.adj)[2]/10 + coef(lm.smb.d15N.adj)[1])-1)) #d15NNO3 at 10 kg m-2 a-1
lm.coef.adj <- c(round(c(lm.adj.sum$coef[2],lm.adj.sum$coef[4],lm.adj.sum$coef[1],lm.adj.sum$coef[3],
                         lm.adj.sum$adj.r.squared), 2),
                 round(c(lm.adj.predict.sum$coef[2],lm.adj.predict.sum$coef[4],lm.adj.predict.sum$coef[1],
                         lm.adj.predict.sum$coef[3], lm.adj.predict.sum$adj.r.squared), 3))

#Only sites with ground observations
ground.obs <- all.nitrate[is.na(all.nitrate$smb.ground) == FALSE,]
lm.smb.d15N.ground <- lm(log(d15N+1)~I(1/smb.ground), data=ground.obs)
lm.smb.predict.d15N.ground <- lm(I(1/smb.ground)~log(d15N+1), data=ground.obs)
lm.ground.sum <- summary(lm.smb.d15N.ground)
lm.ground.predict.sum <- summary(lm.smb.predict.d15N.ground)
smb.range.ground <- c((exp(coef(lm.smb.d15N.ground)[1])-1), #d15NNO3o
                      (exp(coef(lm.smb.d15N.ground)[2]/500 + coef(lm.smb.d15N.ground)[1])-1), #d15NNO3 at 500 kg m-2 a-1
                      (exp(coef(lm.smb.d15N.ground)[2]/200 + coef(lm.smb.d15N.ground)[1])-1), #d15NNO3 at 200 kg m-2 a-1
                      (exp(coef(lm.smb.d15N.ground)[2]/100 + coef(lm.smb.d15N.ground)[1])-1), #d15NNO3 at 100 kg m-2 a-1
                      (exp(coef(lm.smb.d15N.ground)[2]/50 + coef(lm.smb.d15N.ground)[1])-1), #d15NNO3 at 50 kg m-2 a-1
                      (exp(coef(lm.smb.d15N.ground)[2]/20 + coef(lm.smb.d15N.ground)[1])-1), #d15NNO3 at 20 kg m-2 a-1
                      (exp(coef(lm.smb.d15N.ground)[2]/10 + coef(lm.smb.d15N.ground)[1])-1)) #d15NNO3 at 10 kg m-2 a-1
lm.coef.ground <- c(round(c(lm.ground.sum$coef[2],lm.ground.sum$coef[4],lm.ground.sum$coef[1],lm.ground.sum$coef[3],
                            lm.ground.sum$adj.r.squared), 2),
                    round(c(lm.ground.predict.sum$coef[2],lm.ground.predict.sum$coef[4],lm.ground.predict.sum$coef[1],
                            lm.ground.predict.sum$coef[3], lm.ground.predict.sum$adj.r.squared), 3))

#Only sites that lack ground observations
no.ground.obs <- all.nitrate[!is.na(all.nitrate$smb.ground) == FALSE,]
lm.smb.d15N.no.ground <- lm(log(d15N+1)~I(1/smb.mar.adj), data=no.ground.obs)
lm.smb.predict.d15N.no.ground <- lm(I(1/smb.mar.adj)~log(d15N+1), data=no.ground.obs)
lm.no.ground.sum <- summary(lm.smb.d15N.no.ground)
lm.no.ground.predict.sum <- summary(lm.smb.predict.d15N.no.ground)
smb.range.no.ground <- c((exp(coef(lm.smb.d15N.no.ground)[1])-1), #d15NNO3o
                         (exp(coef(lm.smb.d15N.no.ground)[2]/500 + coef(lm.smb.d15N.no.ground)[1])-1), #d15NNO3 at 500 kg m-2 a-1
                         (exp(coef(lm.smb.d15N.no.ground)[2]/200 + coef(lm.smb.d15N.no.ground)[1])-1), #d15NNO3 at 200 kg m-2 a-1
                         (exp(coef(lm.smb.d15N.no.ground)[2]/100 + coef(lm.smb.d15N.no.ground)[1])-1), #d15NNO3 at 100 kg m-2 a-1
                         (exp(coef(lm.smb.d15N.no.ground)[2]/50 + coef(lm.smb.d15N.no.ground)[1])-1), #d15NNO3 at 50 kg m-2 a-1
                         (exp(coef(lm.smb.d15N.no.ground)[2]/20 + coef(lm.smb.d15N.no.ground)[1])-1), #d15NNO3 at 20 kg m-2 a-1
                         (exp(coef(lm.smb.d15N.no.ground)[2]/10 + coef(lm.smb.d15N.no.ground)[1])-1)) #d15NNO3 at 10 kg m-2 a-1
lm.coef.no.ground <- c(round(c(lm.no.ground.sum$coef[2],lm.no.ground.sum$coef[4],lm.no.ground.sum$coef[1],lm.no.ground.sum$coef[3],
                               lm.no.ground.sum$adj.r.squared), 2),
                       round(c(lm.no.ground.predict.sum$coef[2],lm.no.ground.predict.sum$coef[4],lm.no.ground.predict.sum$coef[1],
                               lm.no.ground.predict.sum$coef[3], lm.no.ground.predict.sum$adj.r.squared), 3))

smb.range <- round(data.frame(rbind(smb.range.mix, smb.range.adj, smb.range.ground, smb.range.no.ground)), 1)
colnames(smb.range) <- c("d15NNO3o", "SMB500", "SMB200", "SMB100", "SMB50", "SMB20", "SMB10")
lm.coef <- data.frame(rbind(lm.coef.mix, lm.coef.adj, lm.coef.ground, lm.coef.no.ground))
colnames(lm.coef) <- c("slope", "slope.err", "int", "int.err", "r2", "slope.trns", "slope.trns.err",
                       "int.trns", "int.trns.err", "r2.trns")

#Paternoster et al., 1998 regression parameter comparison
pvalue.regressions <- data.frame(matrix(nrow=2, ncol=8))
rownames(pvalue.regressions) <- c("d15N vs SMB", "SMB vs d15N")
colnames(pvalue.regressions) <- c("slope.mix.adjMAR", "slope.mix.ground", "slope.mix.noground", "slope.ground.noground",
                                  "intercept.mix.adjMAR", "intercept.mix.ground", "intercept.mix.noground", "intercept.ground.noground")

pvalue.regressions[1,1] <- 2*pnorm(-abs((lm.smb.d15N.mix$coef[2]-lm.smb.d15N.ground$coef[2])/sqrt(summary(lm.smb.d15N.mix)$coef[4]+summary(lm.smb.d15N.ground)$coef[4])))
pvalue.regressions[1,2] <- 2*pnorm(-abs((lm.smb.d15N.mix$coef[2]-lm.smb.d15N.no.ground$coef[2])/sqrt(summary(lm.smb.d15N.mix)$coef[4]+summary(lm.smb.d15N.no.ground)$coef[4])))
pvalue.regressions[1,3] <- 2*pnorm(-abs((lm.smb.d15N.mix$coef[2]-lm.smb.d15N.adj$coef[2])/sqrt(summary(lm.smb.d15N.mix)$coef[4]+summary(lm.smb.d15N.adj)$coef[4])))
pvalue.regressions[1,4] <- 2*pnorm(-abs((lm.smb.d15N.no.ground$coef[2]-lm.smb.d15N.ground$coef[2])/sqrt(summary(lm.smb.d15N.no.ground)$coef[4]+summary(lm.smb.d15N.ground)$coef[4])))
pvalue.regressions[1,5] <- 2*pnorm(-abs((lm.smb.d15N.mix$coef[1]-lm.smb.d15N.ground$coef[1])/sqrt(summary(lm.smb.d15N.mix)$coef[3]+summary(lm.smb.d15N.ground)$coef[3])))
pvalue.regressions[1,6] <- 2*pnorm(-abs((lm.smb.d15N.mix$coef[1]-lm.smb.d15N.no.ground$coef[1])/sqrt(summary(lm.smb.d15N.mix)$coef[3]+summary(lm.smb.d15N.no.ground)$coef[3])))
pvalue.regressions[1,7] <- 2*pnorm(-abs((lm.smb.d15N.mix$coef[1]-lm.smb.d15N.adj$coef[1])/sqrt(summary(lm.smb.d15N.mix)$coef[3]+summary(lm.smb.d15N.adj)$coef[3])))
pvalue.regressions[1,8] <- 2*pnorm(-abs((lm.smb.predict.d15N.no.ground$coef[1]-lm.smb.predict.d15N.ground$coef[1])/sqrt(summary(lm.smb.predict.d15N.no.ground)$coef[3]+summary(lm.smb.predict.d15N.ground)$coef[3])))
pvalue.regressions[2,1] <- 2*pnorm(-abs((lm.smb.predict.d15N.mix$coef[2]-lm.smb.predict.d15N.ground$coef[2])/sqrt(summary(lm.smb.predict.d15N.mix)$coef[4]+summary(lm.smb.predict.d15N.ground)$coef[4])))
pvalue.regressions[2,2] <- 2*pnorm(-abs((lm.smb.predict.d15N.mix$coef[2]-lm.smb.predict.d15N.no.ground$coef[2])/sqrt(summary(lm.smb.predict.d15N.mix)$coef[4]+summary(lm.smb.predict.d15N.no.ground)$coef[4])))
pvalue.regressions[2,3] <- 2*pnorm(-abs((lm.smb.predict.d15N.mix$coef[2]-lm.smb.predict.d15N.adj$coef[2])/sqrt(summary(lm.smb.predict.d15N.mix)$coef[4]+summary(lm.smb.predict.d15N.adj)$coef[4])))
pvalue.regressions[2,4] <- 2*pnorm(-abs((lm.smb.predict.d15N.no.ground$coef[2]-lm.smb.predict.d15N.ground$coef[2])/sqrt(summary(lm.smb.predict.d15N.no.ground)$coef[4]+summary(lm.smb.predict.d15N.ground)$coef[4])))
pvalue.regressions[2,5] <- 2*pnorm(-abs((lm.smb.predict.d15N.mix$coef[1]-lm.smb.predict.d15N.ground$coef[1])/sqrt(summary(lm.smb.predict.d15N.mix)$coef[3]+summary(lm.smb.predict.d15N.ground)$coef[3])))
pvalue.regressions[2,6] <- 2*pnorm(-abs((lm.smb.predict.d15N.mix$coef[1]-lm.smb.predict.d15N.no.ground$coef[1])/sqrt(summary(lm.smb.predict.d15N.mix)$coef[3]+summary(lm.smb.predict.d15N.no.ground)$coef[3])))
pvalue.regressions[2,7] <- 2*pnorm(-abs((lm.smb.predict.d15N.mix$coef[1]-lm.smb.predict.d15N.adj$coef[1])/sqrt(summary(lm.smb.predict.d15N.mix)$coef[3]+summary(lm.smb.predict.d15N.adj)$coef[3])))
pvalue.regressions[2,8] <- 2*pnorm(-abs((lm.smb.predict.d15N.no.ground$coef[1]-lm.smb.predict.d15N.ground$coef[1])/sqrt(summary(lm.smb.predict.d15N.no.ground)$coef[3]+summary(lm.smb.predict.d15N.ground)$coef[3])))

#Residuals
lm.smb.d15N.mix.rsd <- log(all.nitrate$d15N+1) - (lm.smb.d15N.mix$coef[2] * 1/all.nitrate$smb.mix + lm.smb.d15N.mix$coef[1])

#Calculating RMSE of function
RMSE <- function(error) {
  sqrt(mean(error^2))
}
RMSE.fit <- lm(I(1/smb.mix)~log(d15N+1), data=all.nitrate)
RMSE(RMSE.fit$residuals)

#####=====END CALCULATING REGRESSIONS #####

#####=====d15N vs SMB PLOTS !!DO NOT RUN UNLESS WANTING PDF CREATION!!=====####

#d15N vs SMB by campaign group
quartz(height=5, width=8)
ggplot(data=all.nitrate) +
  theme_classic() + 
  geom_smooth(aes(x=1/smb.mix, y=log(d15N+1)),
              method='lm', formula=y~x, linetype='longdash', color="gray80", fill="gray80") +
  geom_point(aes(x=1/smb.mix, y=log(d15N+1), color=transect, shape=type), alpha=0.5, size=2) +
  coord_cartesian(xlim = 1/c(22,500), ylim = log(c(-0.020,0.360)+1)) +
  scale_x_reverse(breaks=1/c(25,33,50,100,500),
                  labels=function(x) round(1/x)) + 
  scale_y_continuous(breaks=log(c(0,0.100,0.200,0.300)+1),
                     labels=function(y) round((exp(y)-1)*1000)) +
  scale_color_manual(values = transect.palette) +
  scale_shape_manual(breaks = c("core", "dl", "pit"), values = c(19, 17, 15)) +
  labs(x="SMB (kg m-2 yr-1)", y="d15NNO3arc (‰ vs N2-Air)") +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))

#SMB vs d15N by campaign group
quartz(height=5, width=8)
ggplot(data=all.nitrate) +
  theme_classic() + 
  geom_smooth(aes(y=1/smb.mix, x=log(d15N+1)),
              method='lm', formula=y~x, linetype='longdash', color="gray80", fill="gray80") +
  geom_point(aes(y=1/smb.mix, x=log(d15N+1), color=transect, shape=type), alpha=0.5, size=2) +
  coord_cartesian(xlim = log(c(-20,360)/1000+1), ylim = 1/c(22,500)) +
  scale_y_reverse(breaks=1/c(25,33,50,100,500),
                  labels=function(y) round(1/y)) + 
  scale_x_continuous(breaks=log(c(0,100,200,300)/1000+1),
                     labels=function(x) round(1000*(exp(x)-1))) +
  scale_color_manual(values = transect.palette) +
  scale_shape_manual(values = c(19, 17, 15)) +
  labs(y="SMB (kg m-2 yr-1)", x="d15NNO3arc (‰ vs N2-Air)") +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))

##Outputting PDF Plots of all SMB variants
nitrate.datasets <- list(all.nitrate, all.nitrate, ground.obs, no.ground.obs)
smb.source <- c("smb.mix", "smb.mar.adj", "smb.ground", "smb.mar.adj")

for (i in 1:length(nitrate.datasets)) {
  #d15N vs SMB by campaign group
  pdf(paste("Figures/test/d15N vs SMB by Campaign ", i, ".pdf", sep=""), height=6, width=7.5)
  plot.iter <-  ggplot(data=nitrate.datasets[[i]]) +
    theme_classic() + 
    geom_smooth(aes(x=1/.data[[smb.source[i]]], y=log(d15N+1)),
                method='lm', formula=y~x, linetype='longdash', color="gray65", fill="gray65") +
    geom_point(aes(x=1/.data[[smb.source[i]]], y=log(d15N+1), color=transect, shape=type), alpha=0.5, size=2) +
    coord_cartesian(xlim = 1/c(22,500), ylim = log(c(-20,360)/1000+1)) +
    scale_x_reverse(breaks=1/c(25,33,50,100,500),
                    labels=function(x) round(1/x)) + 
    scale_y_continuous(breaks=log(c(0,100,200,300)/1000+1),
                       labels=function(y) round(1000*(exp(y)-1))) +
    scale_color_manual(values = transect.palette, drop=FALSE) +
    scale_shape_manual(breaks = c("core", "dl", "pit"), values = c(16, 17, 15), drop=FALSE) +
    labs(x="SMB (kg m-2 yr-1)", y="d15NNO3arc (‰ vs N2-Air)") +
    theme(axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(size=16),
          axis.text.x = element_text(size=12),
          plot.title = element_text(hjust=0.5, size=20))
  print(plot.iter)
  dev.off()
  
  
  pdf(paste("Figures/test/SMB vs d15N by Campaign ", i, ".pdf", sep=""), height=6, width=7.5)
  plot.iter <-  ggplot(data=nitrate.datasets[[i]]) +
    theme_classic() + 
    geom_smooth(aes(y=1/.data[[smb.source[i]]], x=log(d15N+1)),
                method='lm', formula=y~x, linetype='longdash', color="gray65", fill="gray65") +
    geom_point(aes(y=1/.data[[smb.source[i]]], x=log(d15N+1), color=transect, shape=type), alpha=0.5, size=2) +
    coord_cartesian(ylim = 1/c(22,500), xlim = log(c(-20,360)/1000+1)) +
    scale_y_reverse(breaks=1/c(25,33,50,100,500),
                    labels=function(y) round(1/y)) + 
    scale_x_continuous(breaks=log(c(0,100,200,300)/1000+1),
                       labels=function(x) round(1000*(exp(x)-1))) +
    scale_color_manual(values = transect.palette, drop=FALSE) +
    scale_shape_manual(breaks = c("core", "dl", "pit"), values = c(16, 17, 15), drop=FALSE) +
    labs(y="SMB (kg m-2 yr-1)", x="d15NNO3arc (‰ vs N2-Air)") +
    theme(axis.title.y = element_text(size=16),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(size=16),
          axis.text.x = element_text(size=12),
          plot.title = element_text(hjust=0.5, size=20))
  print(plot.iter)
  dev.off()
}

#####=====END D15N VS SMB PLOTS #####


##### ABN plotting #####
#==Primary plots and calculations
#d15N reconstructions using inverse function of best mix SMB regression
abn.smb.d15N <- 1/((log(abn$d15N+1)-lm.smb.d15N.mix$coef[1])/lm.smb.d15N.mix$coef[2])

#Mean SMBs per period
abn.smb.20th <- round(c(mean(abn.smb.d15N[abn$yr.btm > 1900], na.rm=T),
                        sd(abn.smb.d15N[abn$yr.btm > 1900], na.rm=T)), 0) #20th century
abn.smb.pre20th <- round(c(mean(abn.smb.d15N[abn$yr.top <= 1900], na.rm=T),
                           sd(abn.smb.d15N[abn$yr.top <= 1900], na.rm=T)), 0) #Pre-20th century
abn.smb.17th <- round(c(mean(abn.smb.d15N[abn$yr.top < 1700 & abn$yr.btm >= 1600], na.rm=T),
                        sd(abn.smb.d15N[abn$yr.top < 1700 & abn$yr.btm >= 1600], na.rm=T)), 0) #1600s CE

#Coefficient of variation (sd/mean)
sd(abn.smb.d15N, na.rm=TRUE)/mean(abn.smb.d15N, na.rm=TRUE)

#Matching NO3 resolution of ABN to density record
abn.density.iter <- NA
abn.smb.iter <- NA
gpr.upstream.iter <- NA
for(i in 2:nrow(abn.density)+1) {
  if(length(abn.smb.d15N[abn$yr.top.BP <= abn.density$yr.top.BP[i] &
                         abn$yr.top.BP > abn.density$yr.top.BP[i-1]]) > 0) {
    abn.smb.iter[i-1] <- mean(abn.smb.d15N[abn$yr.top.BP <= abn.density$yr.top.BP[i] &
                                             abn$yr.top.BP > abn.density$yr.top.BP[i-1]], na.rm=TRUE)
    gpr.upstream.iter[i-1] <- mean(abn$smb.gpr.upstream[abn$yr.top.BP <= abn.density$yr.top.BP[i] &
                                                          abn$yr.top.BP > abn.density$yr.top.BP[i-1]], na.rm=TRUE)
  } else {
    abn.smb.iter[i-1] <- NA
    gpr.upstream.iter[i-1] <- NA
  }
}

#Detrending SMBd15N based on SMBGPR
abn.detrend <- abn.smb.d15N - abn$smb.gpr.upstream
gpr.offset <- mean(abn$smb.gpr.upstream, na.rm=TRUE) - mean(abn.smb.d15N, na.rm=TRUE) #calculating the SMB-offset between GPR and d15N estimates

#Table of correlations for SMB reconstructions
abn.smb.cor <- data.frame(matrix(ncol=2, nrow=3)) #Pearson
rownames(abn.smb.cor) <- c("d15N", "Density", "GPR")
colnames(abn.smb.cor) <- c("Density", "GPR")
abn.smb.cor.spearman <- data.frame(matrix(ncol=2, nrow=3)) #Spearman
rownames(abn.smb.cor.spearman) <- c("d15N", "Density", "GPR")
colnames(abn.smb.cor.spearman) <- c("Density", "GPR")
abn.smb.cor[1,1] <- cor.test(abn.density$smb, abn.smb.iter, use="complete.obs")$estimate
abn.smb.cor.spearman[1,1] <- cor.test(abn.density$smb, abn.smb.iter, method="spearman", use="complete.obs")$estimate
abn.smb.cor[1,2] <- cor.test(gpr.upstream.iter, abn.smb.iter, use="complete.obs")$estimate
abn.smb.cor.spearman[1,2] <- cor.test(gpr.upstream.iter, abn.smb.iter, method="spearman", use="complete.obs")$estimate
abn.smb.cor[2,2] <- cor.test(abn.density$smb, gpr.upstream.iter, use="complete.obs")$estimate
abn.smb.cor.spearman[2,2] <- cor.test(abn.density$smb, gpr.upstream.iter, method="spearman", use="complete.obs")$estimate
round(abn.smb.cor, 2)
round(abn.smb.cor.spearman, 2)

#Calculating 50 yr running means and sds
yr.iter.50 <- NA
abn.smb.d15N.50yr <- NA
abn.detrend.50yr <- NA
abn.smb.d15Nsd.50yr <- NA
abn.detrendsd.50yr <- NA
for (i in 1:nrow(abn)) {
  yr.iter.50[i] <- abn$yr.top.BP[i]
  if (yr.iter.50[i] < 600) {
    abn.smb.d15N.50yr[i] <- mean(abn.smb.d15N[abn$yr.top.BP > yr.iter.50[i] & abn$yr.top.BP <= (yr.iter.50[i]+49)], na.rm = TRUE)
    abn.detrend.50yr[i] <- mean(abn.detrend[abn$yr.top.BP > yr.iter.50[i] & abn$yr.top.BP <= (yr.iter.50[i]+49)], na.rm = TRUE)
    abn.smb.d15Nsd.50yr[i] <- sd(abn.smb.d15N[abn$yr.top.BP > yr.iter.50[i] & abn$yr.top.BP <= (yr.iter.50[i]+49)], na.rm = TRUE)
    abn.detrendsd.50yr[i] <- sd(abn.detrend[abn$yr.top.BP > yr.iter.50[i] & abn$yr.top.BP <= (yr.iter.50[i]+49)], na.rm = TRUE)
  }
  else {
    abn.smb.d15N.50yr[i] <- NA    
    abn.detrend.50yr[i] <- NA
    abn.smb.d15Nsd.50yr[i] <- NA    
    abn.detrendsd.50yr[i] <- NA
  }
}

#=Plots
#Plotting SMB prediction from d15N
p.abn.smb.d15N <- ggplot() +
  theme_classic() + 
  geom_ribbon(aes(x=yr.iter.50+25, ymax=abn.smb.d15N.50yr+abn.smb.d15Nsd.50yr, ymin=abn.smb.d15N.50yr-abn.smb.d15Nsd.50yr),
              fill=adjustcolor("firebrick1", alpha=0.4)) +
  geom_step(aes(x=abn$yr.top.BP, y=abn.smb.d15N), color="coral1") +
  geom_line(aes(x=yr.iter.50+25, y=abn.smb.d15N.50yr), color="firebrick", lwd=1.5) +
  scale_x_continuous() +
  scale_y_continuous(sec.axis = sec_axis(trans=~exp(.^-1*lm.smb.d15N.mix$coef[2]+lm.smb.d15N.mix$coef[1])-1,
                                         breaks=c(0.025,0.05,0.1,0.15), labels=c(25, 50, 100, 150))) +
  labs(x="Year BP", y="SMB (kg m-2 yr-1)", 
       title="ABN: Predicted SMB from d15N") +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))
#pdf("Figures/abn_smb_d15N.pdf", height=5, width=8)
quartz(height=5, width=8)
p.abn.smb.d15N
#dev.off()

#Comparing d15N to density SMB reconstructions, resolution matched 
p.abn.smb.compare <- ggplot() +
  theme_classic() + 
  geom_step(aes(x=abn.density$yr.top.BP, y=abn.smb.iter), color="coral1") +
  geom_step(aes(x=abn.density$yr.top.BP, y=abn.density$smb), color="gray60") +
  geom_smooth(aes(x=abn.density$yr.top.BP, y=abn.density$smb), method="loess", color="gray40",
              fill="gray40", span=0.33, alpha=0.2) +
  geom_smooth(aes(x=abn.density$yr.top.BP, y=abn.smb.iter), method="loess", color="firebrick",
              fill="firebrick", span=0.33, alpha=0.2) +
  geom_line(aes(x=abn$yr.top.BP, y=abn$smb.gpr.upstream), color="dodgerblue4", lwd = 1, lty=2) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(x="Year BP", y="SMB (kg m-2 yr-1)", 
       title="ABN: Predicted SMB from d15N, Density, and GPR") +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))
#pdf("Figures/abn_smb_compare.pdf", height=5, width=8)
quartz(height=5, width=8)
p.abn.smb.compare
#dev.off()

#=ABN SMBs Detrended based on GPR upstream observations
p.abn.smb.detrend <- ggplot() +
  theme_classic() + 
  geom_hline(yintercept=0, linetype="dashed") +
  geom_ribbon(aes(x=yr.iter.50+25, ymax=abn.detrend.50yr+abn.detrendsd.50yr+gpr.offset, ymin=abn.detrend.50yr-abn.detrendsd.50yr+gpr.offset),
              fill=adjustcolor("firebrick1", alpha=0.4)) +
  geom_step(aes(x=abn$yr.top.BP, y=abn.detrend+gpr.offset), color="coral1") +
  geom_line(aes(x=yr.iter.50+25, y=abn.detrend.50yr+gpr.offset), color="firebrick", lwd=1.5) +
  scale_x_continuous() +
  scale_y_continuous() +
  labs(x="Year BP", y="SMB anomaly (kg m-2 a-1)", 
       title="ABN: Predicted SMB from d15N, Upstream GPR Detrended") +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))
#pdf("Figures/abn_smb_gpr_detrend.pdf", height=5, width=8)
quartz(height=5, width=8)
p.abn.smb.detrend
#dev.off()

#Three plots combined
plot.list <- list(p.abn.smb.d15N, p.abn.smb.compare, p.abn.smb.detrend)
#pdf("Figures/abn_smb_3plot.pdf", height=10, width=8)
quartz(height=10, width=8)
plot_grid(plotlist=plot.list, ncol=1)
#dev.off()

#Determining changes to SMBd15N vs SMBdensity correlations with aggregration
yr.iter <- list()
abn.smb.d15N.iter <- list()
abn.smb.density.iter <- list()
abn.smb.d15Nsd.iter <- list()
cor.iter.abn <- NA
#Calculating multiyear running means and sds
for(j in 1:200) { #j is the number of years aggregated
  yr.iter[[j]] <- NA
  abn.smb.d15N.iter[[j]] <- NA
  abn.smb.density.iter[[j]] <- NA
  abn.smb.d15Nsd.iter[[j]] <- NA
  for (i in 1:nrow(abn)) {
    yr.iter[[j]][i] <- abn$yr.top.BP[i]
    if (yr.iter[[j]][i] < max(abn$yr.top.BP-j)) { #Aggregating by the year interval set by j
      abn.smb.d15N.iter[[j]][i] <- mean(abn.smb.iter[abn.density$yr.top.BP > yr.iter[[j]][i] & abn.density$yr.top.BP <= (yr.iter[[j]][i]+(j))], na.rm = TRUE)
      abn.smb.density.iter[[j]][i] <- mean(abn.density$smb[abn.density$yr.top.BP > yr.iter[[j]][i] & abn.density$yr.top.BP <= (yr.iter[[j]][i]+(j))], na.rm = TRUE)
    }
    else {
      abn.smb.d15N.iter[[j]][i] <- NA    
      abn.smb.density.iter[[j]][i] <- NA    
    }
  }
  cor.iter.abn[j] <-cor(abn.smb.d15N.iter[[j]],abn.smb.density.iter[[j]],use="complete.obs")
}
abn.smb.cor.yr.iter <- data.frame(seq(1,200),cor.iter.abn)
colnames(abn.smb.cor.yr.iter) <- c("years.aggregated", "r")
print(abn.smb.cor.yr.iter) #Note that the 1 and 2 yr aggregations are lower resolution than the actual sample resolution and should be ignored
quartz()
plot(seq(1,200),cor.iter.abn, type="l", main="ABN Correlation with Aggregation of Years")

##### END ABN #####


#####=========================MAR ADJUSTMENT CALCULATIONS============================#####
# This script corrects the biases in MAR ERA SMB estimates (Agosta 2019) for EAA based on comparison
#   between the MAR estimates and actual direct ground measurements.
# Written by Pete D Akers, Dec 2020

mar.orig <- read.csv("mar_adjust_mar.csv", sep=";", header=TRUE)
mar.direct.all <- read.csv("mar_adjust_direct_all.csv", sep=";", header=TRUE)
mar.direct.small <- mar.direct.all[mar.direct.all$less175 == 1,] #All sites less than 175 in either MAR or direct
mar.direct.big <- mar.direct.all[mar.direct.all$more110 == 1,] #All sites greater than 110 in MAR or direct

lm.mar.adjust.small <- lm(smb.direct~smb.mar.orig, data=mar.direct.small)
lm.mar.adjust.big <- lm(smb.direct~smb.mar.orig, data=mar.direct.big)

#Find the intersection of the regressions
cm <- rbind(coef(lm.mar.adjust.small),coef(lm.mar.adjust.big)) # Coefficient matrix
c(-solve(cbind(cm[,2],-1)) %*% cm[,1])

#The regressions cross at x=138, so use that as split
mar.orig.small <- mar.orig[mar.orig$smb.mar.orig <= 138,]
mar.orig.big <- mar.orig[mar.orig$smb.mar.orig > 138,] #all nonculled sites with smb greater than 175
predict.small <- cbind(mar.orig.small, data.frame(predict.lm(lm.mar.adjust.175, newdata=mar.orig.small, interval="prediction")))
predict.big <- cbind(mar.orig.big, data.frame(predict.lm(lm.mar.adjust.big, newdata=mar.orig.big, interval="prediction")))
predict.all <- rbind(predict.small, predict.big)

write.csv(predict.all, "mar_adjust_results.csv", row.names=FALSE)

#Split Plot
#pdf("Figures/MAR vs Direct SMB.pdf", height=6, width=6)
quartz()
ggplot() +
  theme_classic() + 
  geom_smooth(aes(x=mar.direct.175$smb.mar.orig, y=mar.direct.175$smb.direct), method="lm", color="firebrick",
              fill=adjustcolor("firebrick4", alpha=0.2)) +
  geom_smooth(aes(x=mar.direct.big$smb.mar.orig, y=mar.direct.big$smb.direct), method="lm", color="dodgerblue4",
              fill=adjustcolor("dodgerblue4", alpha=0.2)) +  
  geom_point(aes(x=mar.direct.175$smb.mar.orig, y=mar.direct.175$smb.direct), col="firebrick", shape=16, alpha=0.6) +
  geom_point(aes(x=mar.direct.big$smb.mar.orig, y=mar.direct.big$smb.direct), col="dodgerblue4", shape=5, alpha=0.6) +
  geom_abline(slope=1,intercept=0) +  
  geom_vline(xintercept=139, linetype=3) +  
  labs(x="MAR SMB (kg/m2/yr)", y="Observed SMB (kg/m2/yr)") +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))
#dev.off()

##### END MAR ADJUSTMENT #####

#####=======DomeA and WAIS Reconstructions==========#####
da <- read.csv("DA_iso_ion.csv", header=TRUE, sep=";") #Jiang et al., 2019: DA2005 ice core nitrate
da.smb <- read.csv("DA_smb.csv", header=TRUE, sep=";") #Jiang et al., 2012: DA2005 ice core density
da.smb.d15N <- 1/((log(da$d15N+1)-lm.smb.d15N.mix$coef[1])/lm.smb.d15N.mix$coef[2]) #DA2005 SMBd15N reconstruction
wais <- read.csv("wais_d15N.csv", header=TRUE, sep=";") #Sofen et al., 2014: WD06A ice core nitrate and density
wais.smb.d15N <- 1/((log(wais$d15N+1)-lm.smb.d15N.mix$coef[1])/lm.smb.d15N.mix$coef[2]) #WD06A SMBd15N reconstruction
wais.smb.d15N.ground <- 1/((log(wais$d15N+1)-lm.smb.d15N.ground$coef[1])/lm.smb.d15N.ground$coef[2]) #WD06A SMBd15N reconstruction based on ground observations only regression

#Matching density resolution of DA to d15N record
da.smb.iter <- NA
for(i in 1:length(da.smb.d15N)) {
  if(length(da.smb$smb[da.smb$yr.top.BP <= da$yr.top.BP[i+1] &
                       da.smb$yr.top.BP > da$yr.top.BP[i]]) > 0) {
    da.smb.iter[i] <- mean(da.smb$smb[da.smb$yr.top.BP <= da$yr.top.BP[i+1] &
                                        da.smb$yr.top.BP > da$yr.top.BP[i]], na.rm=TRUE)
  } else {
    da.smb.iter[i] <- NA
  }
  if(i==length(da.smb.d15N)){
    da.smb.iter[i] <- mean(da.smb$smb[da.smb$yr.top.BP <= da$yr.btm.BP[i] &
                                        da.smb$yr.top.BP > da$yr.top.BP[i]], na.rm=TRUE)
  }
}

#DA2005 SMBd15N reconstruction
#pdf("Figures/DA SMB d15N Compare.pdf", height=5, width=8)
quartz(height=5, width=8)
ggplot()+
  theme_classic()+
  geom_step(aes(x=da$yr.top.BP, y=da.smb.d15N), color="darkgreen") +
  geom_step(aes(x=da$yr.top.BP, y=da.smb.iter), color="gray60") +
  labs(x="Year BP", y="SMB (kg m-2 a-1)") +
  scale_y_continuous(sec.axis = sec_axis(trans=~exp(.^-1*lm.smb.d15N.mix$coef[2]+lm.smb.d15N.mix$coef[1])-1,
                                         breaks=c(0.25,0.275,0.3,0.325), labels=c(250,275,300,325))) +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))
#dev.off()

#WD06A SMBd15N reconstruction
#pdf("Figures/WAIS SMB d15N Compare.pdf", height=5, width=8)
quartz(height=5, width=8)
ggplot()+
  theme_classic()+
  geom_step(aes(x=wais$yr.top.BP, y=wais.smb.d15N), color="mediumorchid4") +
  geom_step(aes(x=wais$yr.top.BP, y=wais.smb.d15N.ground), color="hotpink3") +
  geom_step(aes(x=wais$yr.top.BP, y=wais$smb), color="gray60") +
  labs(x="Year BP", y="SMB (kg m-2 a-1)") +
  scale_y_continuous(sec.axis = sec_axis(trans=~exp(.^-1*lm.smb.d15N.mix$coef[2]+lm.smb.d15N.mix$coef[1])-1,
                                         breaks=c(0.002,0.004,0.006,0.008,0.010,0.012, 0.014, 0.016), labels=c(2, 4, 6, 8, 10, 12, 14, 16))) +
  theme(axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20))
#dev.off()
