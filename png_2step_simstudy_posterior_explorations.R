library(INLA); library(raster); library(maptools)
library(gtools); library(sp); library(spdep); library(rgdal)
library(fields); library(mvtnorm);  library(geoR)
library(actuar);library(viridisLite);require(grid);require(gridExtra)
require(lattice);require(tidyverse);require(MASS);library(tmap)
library(tmaptools);library(sf)
#install.packages("cartography")

library(cartography) # mapping dedicated package
#install.packages("OpenStreetMap")
library(OpenStreetMap)

#This script if for 9 areas

path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/Working/Chris/reviews/sim_study2"
out_path <- paste0(path, "/outputs")
file_path <- paste0(path, "/file")
setwd(out_path)

pop.cover <- c(100, 80, 60, 40, 20)
bld.cover <- c(100, 95, 90, 85, 80, 75, 70, 65)
metric <- c("mae", "rmse", "bias", "corr")
method <- c("onestep", "twostep")

n.pop.cover <- length(pop.cover)
n.bld.cover <- length(bld.cover)
n.metric <- length(metric)
n.method <- length(method)

##---build the dataframe for the metrics
n.size <- n.pop.cover*n.bld.cover*n.metric*n.method
dim(dat.met <- data.frame(expand.grid(method=method,
                                      bld_cover=bld.cover,
                                      pop_cover=pop.cover,
                                      metric=metric)))

###

#paste0(out_path, "/outputs_for_80%_pop_count/95%_bldg_count")

#dat_met <- matrix(0, nrow=n.size, ncol=5)
dat_met1 <- list()
dat_met2 <- list()
for(j in 1:n.pop.cover)
{
  pathp <- paste0(out_path,"/outputs_for_", pop.cover[j],"%","_pop_count")
  for(k in 1:n.bld.cover)
  {
    pathb <- paste0(pathp,"/", bld.cover[k], "%","_bldg_count")
    met0 <- read.csv(paste0(pathb, "/fit_metrics_0.csv"))
    met1 <- read.csv(paste0(pathb, "/fit_metrics_1.csv"))
    met0[1] <- bld.cover[k]
    met1[1] <- bld.cover[k]
    met0 <- c(pop.cover[j], met0)
    met1 <- c(pop.cover[j], met1)
    dat_met1[[k]] = rbind(met0, met1) 
  }
  dat_met2[[j]] = dat_met1
} 
#dat_met2
unnest_list <- unlist(dat_met2, recursive = FALSE)  #--unnest the list
dim(metrics <- as.data.frame(matrix(unlist(do.call(rbind, unnest_list)),
                                    nrow=80, ncol=6)))
##
names(metrics) <- c("pop_cover", "bld_cover", "mae","rmse", "bias", "corr") #-rename columns
metrics$method <- rep(c("onestep", "twosteps"),nrow(metrics)/2)#--add 'method' col

head(metrics)
write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE)

# Convert to long format for plotting 
require(reshape2)
require(ggpubr)
dim(met_long <- melt(metrics, id.vars=c("pop_cover","bld_cover", "method"),
                     value.name="estimate", variable.name = "metric"))

met_long$method = factor(met_long$method)
met_long$pop_cover = factor(met_long$pop_cover)
head(met_long)
table(met_long$metric)
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

##---Plots
variable_names <- list(
  "mae" = "Mean  Absolute \n Error (MAE)" ,
  "rmse" = "Root Mean Square \n Error (RMSE)",
  "bias" = "Absolute Bias",
  "corr" = "Correlation \n Coefficient"
)

levels(met_long$method) <- c("BHM","TSBHM")
grp <- levels(met_long$method)

variable_labeller2 <- function(variable,value){
  if (variable=='metric') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}



plot_metrics <- ggplot(met_long, aes(x=bld_cover, y=estimate))+
  geom_point(aes(colour=pop_cover, shape=pop_cover), size=2)+
  geom_line(aes(colour=pop_cover), size=1)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_grid(metric~method, scales="free", 
             space="free_x",  labeller= variable_labeller2)
plot_metrics

plot_met <- ggpar(plot_metrics, xlab="Sattelite Observation Coverage (%)", ylab="Estimate",
               legend = "right", legend.title = "Survey \n Coverage (%)",
               palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0D0887FF", "#993333"),
               #colour = "bld_cover",
               #shape= "bld_cover",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(16),
               font.y = c(16),
               xtickslab.rt = 45, ytickslab.rt = 45)

plot_met
##


####============------------------------
#---------------Make Scatter plots for pop counts------------
dat_cu0<- list()
dat_cu1 <- list()
for(j in 1:n.pop.cover)
{
  pathp <- paste0(out_path,"/outputs_for_", pop.cover[j],"%","_pop_count")
  for(k in 1:n.bld.cover)
  {
    pathb <- paste0(pathp,"/", bld.cover[k], "%","_bldg_count")
    cu0 <- read.csv(paste0(pathb, "/CU_estimates_0.csv"))
    cu1 <- read.csv(paste0(pathb, "/CU_estimates_1.csv"))
    names(cu1); names(cu0)
    cu0$mean1 <- cu1$mean
    cu0$lower1 <- cu1$lower
    cu0$upper1 <- cu1$upper
    cu0$sd1 <- cu1$sd_pop_hat
    cu0$dens1 <- cu1$mean_dens_hat
    cu0$upper1 <- cu1$upper
    cu0$sd1 <- cu1$sd_pop_hat
    cu0$bld_cover <- rep(bld.cover[k], nrow(cu0))
    cu0$pop_cover <- rep(pop.cover[j], nrow(cu0))
    #cu0$method <- rep(c("onestep", "twosteps"),nrow(cu0)/2)
    dat_cu0[[k]] = cu0
  }
  dat_cu1[[j]] = dat_cu0
} 
#dat_cu1  #- very huge and 

##--Extract 
datpop100 <- dat_cu1[[1]]
datpop80 <- dat_cu1[[2]]
datpop60 <- dat_cu1[[3]]
datpop40 <- dat_cu1[[4]]
datpop20 <- dat_cu1[[5]]

class(datpop100)

##---Select variables to include
Var2Include <- c("lon", "lat", "prov2_ID", "set_typ2_ID",
                 "x2", "x3", "x4", "x5", "x6", "bld", "pop",
                 "dens", "popm", "bldm", "prov_ranef","mean_dens_hat", "dens1",
                 "mean", "lower", "upper", "mean1", "lower1", "upper1",
                 "pop_cover", "bld_cover")

##########
##----Extract data for 60% survey coverage
dim(p100 <- rbind(datpop100[[1]][,Var2Include],
                  datpop100[[2]][,Var2Include],
                  datpop100[[3]][,Var2Include],
                  datpop100[[4]][,Var2Include],
                  datpop100[[5]][,Var2Include],
                  datpop100[[6]][,Var2Include],
                  datpop100[[7]][,Var2Include],
                  datpop100[[8]][,Var2Include]))



dim(p100_0 <- p100[,names(p100)[c(1:16,18:20,24,25)]])
p100_0$method <- rep("BHM", nrow(p100_0))
p100_0$method  <- factor(p100_0$method)
p100_0$bld_cover  <- factor(p100_0$bld_cover)

####
dim(p100_1 <- p100[,names(p100)[c(1:15,17,21:25)]])
p100_1$method <- rep("TSBHM", nrow(p100_1))
p100_1$method  <- factor(p100_1$method)
p100_1$bld_cover  <- factor(p100_1$bld_cover)

names(p100_1)= names(p100_0)



###------
dim(p100_0_b1 <- p100_0[p100_0$bld_cover=="100",])
dim(p100_0_b2 <- p100_0[p100_0$bld_cover=="95",])
dim(p100_0_b3 <- p100_0[p100_0$bld_cover=="90",])
dim(p100_0_b4 <- p100_0[p100_0$bld_cover=="85",])
dim(p100_0_b5 <- p100_0[p100_0$bld_cover=="80",])
dim(p100_0_b6 <- p100_0[p100_0$bld_cover=="75",])
dim(p100_0_b7 <- p100_0[p100_0$bld_cover=="70",])
dim(p100_0_b8 <- p100_0[p100_0$bld_cover=="65",])

##
dim(p100_1_b1 <- p100_1[p100_1$bld_cover=="100",])
dim(p100_1_b2 <- p100_1[p100_1$bld_cover=="95",])
dim(p100_1_b3 <- p100_1[p100_1$bld_cover=="90",])
dim(p100_1_b4 <- p100_1[p100_1$bld_cover=="85",])
dim(p100_1_b5 <- p100_1[p100_1$bld_cover=="80",])
dim(p100_1_b6 <- p100_1[p100_1$bld_cover=="75",])
dim(p100_1_b7 <- p100_1[p100_1$bld_cover=="70",])
dim(p100_1_b8 <- p100_1[p100_1$bld_cover=="65",])


#---80%---------------------------------------------------------
##----Extract data for 60% survey coverage
dim(p80 <- rbind(datpop80[[1]][,Var2Include],
                 datpop80[[2]][,Var2Include],
                 datpop80[[3]][,Var2Include],
                 datpop80[[4]][,Var2Include],
                 datpop80[[5]][,Var2Include],
                 datpop80[[6]][,Var2Include],
                 datpop80[[7]][,Var2Include],
                 datpop80[[8]][,Var2Include]))



dim(p80_0 <- p80[,names(p80)[c(1:16,18:20,24,25)]])
p80_0$method <- rep("BHM", nrow(p80_0))
p80_0$method  <- factor(p80_0$method)
p80_0$bld_cover  <- factor(p80_0$bld_cover)

####
dim(p80_1 <- p80[,names(p80)[c(1:15,17,21:25)]])
p80_1$method <- rep("TSBHM", nrow(p80_1))
p80_1$method  <- factor(p80_1$method)
p80_1$bld_cover  <- factor(p80_1$bld_cover)

names(p80_1)= names(p80_0)



###------
dim(p80_0_b1 <- p80_0[p80_0$bld_cover=="100",])
dim(p80_0_b2 <- p80_0[p80_0$bld_cover=="95",])
dim(p80_0_b3 <- p80_0[p80_0$bld_cover=="90",])
dim(p80_0_b4 <- p80_0[p80_0$bld_cover=="85",])
dim(p80_0_b5 <- p80_0[p80_0$bld_cover=="80",])
dim(p80_0_b6 <- p80_0[p80_0$bld_cover=="75",])
dim(p80_0_b7 <- p80_0[p80_0$bld_cover=="70",])
dim(p80_0_b8 <- p80_0[p80_0$bld_cover=="65",])

##
dim(p80_1_b1 <- p80_1[p80_1$bld_cover=="100",])
dim(p80_1_b2 <- p80_1[p80_1$bld_cover=="95",])
dim(p80_1_b3 <- p80_1[p80_1$bld_cover=="90",])
dim(p80_1_b4 <- p80_1[p80_1$bld_cover=="85",])
dim(p80_1_b5 <- p80_1[p80_1$bld_cover=="80",])
dim(p80_1_b6 <- p80_1[p80_1$bld_cover=="75",])
dim(p80_1_b7 <- p80_1[p80_1$bld_cover=="70",])
dim(p80_1_b8 <- p80_1[p80_1$bld_cover=="65",])

###---60%---------------------------------------------------------
##----Extract data for 60% survey coverage
dim(p60 <- rbind(datpop60[[1]][,Var2Include],
                 datpop60[[2]][,Var2Include],
                 datpop60[[3]][,Var2Include],
                 datpop60[[4]][,Var2Include],
                 datpop60[[5]][,Var2Include],
                 datpop60[[6]][,Var2Include],
                 datpop60[[7]][,Var2Include],
                 datpop60[[8]][,Var2Include]))



dim(p60_0 <- p60[,names(p60)[c(1:16,18:20,24,25)]])
p60_0$method <- rep("BHM", nrow(p60_0))
p60_0$method  <- factor(p60_0$method)
p60_0$bld_cover  <- factor(p60_0$bld_cover)

####
dim(p60_1 <- p60[,names(p60)[c(1:15,17,21:25)]])
p60_1$method <- rep("TSBHM", nrow(p60_1))
p60_1$method  <- factor(p60_1$method)
p60_1$bld_cover  <- factor(p60_1$bld_cover)

names(p60_1)= names(p60_0)



###------
dim(p60_0_b1 <- p60_0[p60_0$bld_cover=="100",])
dim(p60_0_b2 <- p60_0[p60_0$bld_cover=="95",])
dim(p60_0_b3 <- p60_0[p60_0$bld_cover=="90",])
dim(p60_0_b4 <- p60_0[p60_0$bld_cover=="85",])
dim(p60_0_b5 <- p60_0[p60_0$bld_cover=="80",])
dim(p60_0_b6 <- p60_0[p60_0$bld_cover=="75",])
dim(p60_0_b7 <- p60_0[p60_0$bld_cover=="70",])
dim(p60_0_b8 <- p60_0[p60_0$bld_cover=="65",])

##
dim(p60_1_b1 <- p60_1[p60_1$bld_cover=="100",])
dim(p60_1_b2 <- p60_1[p60_1$bld_cover=="95",])
dim(p60_1_b3 <- p60_1[p60_1$bld_cover=="90",])
dim(p60_1_b4 <- p60_1[p60_1$bld_cover=="85",])
dim(p60_1_b5 <- p60_1[p60_1$bld_cover=="80",])
dim(p60_1_b6 <- p60_1[p60_1$bld_cover=="75",])
dim(p60_1_b7 <- p60_1[p60_1$bld_cover=="70",])
dim(p60_1_b8 <- p60_1[p60_1$bld_cover=="65",])


###---40%--------------------------------------------------------------
##----Extract data for 40% survey coverage
dim(p40 <- rbind(datpop40[[1]][,Var2Include],
                 datpop40[[2]][,Var2Include],
                 datpop40[[3]][,Var2Include],
                 datpop40[[4]][,Var2Include],
                 datpop40[[5]][,Var2Include],
                 datpop40[[6]][,Var2Include],
                 datpop40[[7]][,Var2Include],
                 datpop40[[8]][,Var2Include]))



dim(p40_0 <- p40[,names(p40)[c(1:16,18:20,24,25)]])
p40_0$method <- rep("BHM", nrow(p40_0))
p40_0$method  <- factor(p40_0$method)
p40_0$bld_cover  <- factor(p40_0$bld_cover)

####
dim(p40_1 <- p40[,names(p40)[c(1:15,17,21:25)]])
p40_1$method <- rep("TSBHM", nrow(p40_1))
p40_1$method  <- factor(p40_1$method)
p40_1$bld_cover  <- factor(p40_1$bld_cover)

names(p40_1)= names(p40_0)



###------
dim(p40_0_b1 <- p40_0[p40_0$bld_cover=="100",])
dim(p40_0_b2 <- p40_0[p40_0$bld_cover=="95",])
dim(p40_0_b3 <- p40_0[p40_0$bld_cover=="90",])
dim(p40_0_b4 <- p40_0[p40_0$bld_cover=="85",])
dim(p40_0_b5 <- p40_0[p40_0$bld_cover=="80",])
dim(p40_0_b6 <- p40_0[p40_0$bld_cover=="75",])
dim(p40_0_b7 <- p40_0[p40_0$bld_cover=="70",])
dim(p40_0_b8 <- p40_0[p40_0$bld_cover=="65",])

##
dim(p40_1_b1 <- p40_1[p40_1$bld_cover=="100",])
dim(p40_1_b2 <- p40_1[p40_1$bld_cover=="95",])
dim(p40_1_b3 <- p40_1[p40_1$bld_cover=="90",])
dim(p40_1_b4 <- p40_1[p40_1$bld_cover=="85",])
dim(p40_1_b5 <- p40_1[p40_1$bld_cover=="80",])
dim(p40_1_b6 <- p40_1[p40_1$bld_cover=="75",])
dim(p40_1_b7 <- p40_1[p40_1$bld_cover=="70",])
dim(p40_1_b8 <- p40_1[p40_1$bld_cover=="65",])


##---20%----------------------------------------------------------------
##----Extract data for 20% survey coverage
dim(p20 <- rbind(datpop20[[1]][,Var2Include],
                 datpop20[[2]][,Var2Include],
                 datpop20[[3]][,Var2Include],
                 datpop20[[4]][,Var2Include],
                 datpop20[[5]][,Var2Include],
                 datpop20[[6]][,Var2Include],
                 datpop20[[7]][,Var2Include],
                 datpop20[[8]][,Var2Include]))


####---BHM
dim(p20_0 <- p20[,names(p20)[c(1:16,18:20,24,25)]])
p20_0$method <- rep("BHM", nrow(p20_0))
p20_0$method  <- factor(p20_0$method)
p20_0$bld_cover  <- factor(p20_0$bld_cover)


##----TSBHM
####
dim(p20_1 <- p20[,names(p20)[c(1:15,17,21:25)]])
p20_1$method <- rep("TSBHM", nrow(p20_1))
p20_1$method  <- factor(p20_1$method)
p20_1$bld_cover  <- factor(p20_1$bld_cover)


names(p20_1)= names(p20_0)



###------
dim(p20_0_b1 <- p20_0[p20_0$bld_cover=="100",])
dim(p20_0_b2 <- p20_0[p20_0$bld_cover=="95",])
dim(p20_0_b3 <- p20_0[p20_0$bld_cover=="90",])
dim(p20_0_b4 <- p20_0[p20_0$bld_cover=="85",])
dim(p20_0_b5 <- p20_0[p20_0$bld_cover=="80",])
dim(p20_0_b6 <- p20_0[p20_0$bld_cover=="75",])
dim(p20_0_b7 <- p20_0[p20_0$bld_cover=="70",])
dim(p20_0_b8 <- p20_0[p20_0$bld_cover=="65",])

##
dim(p20_1_b1 <- p20_1[p20_1$bld_cover=="100",])
dim(p20_1_b2 <- p20_1[p20_1$bld_cover=="95",])
dim(p20_1_b3 <- p20_1[p20_1$bld_cover=="90",])
dim(p20_1_b4 <- p20_1[p20_1$bld_cover=="85",])
dim(p20_1_b5 <- p20_1[p20_1$bld_cover=="80",])
dim(p20_1_b6 <- p20_1[p20_1$bld_cover=="75",])
dim(p20_1_b7 <- p20_1[p20_1$bld_cover=="70",])
dim(p20_1_b8 <- p20_1[p20_1$bld_cover=="65",])

##------------------------------------------------------------------------------------------------
##-----SCATTER PLOTS---------------------------------------------------------------
#################################----------------------------------------------
p100_2 <- p100_0
p100_2$mean <- p100_0$mean
p100_2$method <- rep("FullData", nrow(p100_2))
dim(p100m <- rbind(p100_0, p100_1, p100_2))



p80_2 <- p80_0
p80_2$mean <- p100_0$mean
p80_2$method <- rep("FullData", nrow(p80_2))
dim(p80m <- rbind(p80_0, p80_1, p80_2))



p60_2 <- p60_0
p60_2$mean <- p100_0$mean
p60_2$method <- rep("FullData", nrow(p60_2))
dim(p60m <- rbind(p60_0, p60_1, p60_2))



p40_2 <- p40_0
p40_2$mean <- p100_0$mean
p40_2$method <- rep("FullData", nrow(p40_2))
dim(p40m <- rbind(p40_0, p40_1, p40_2))



p20_2 <- p20_0
p20_2$mean <- p100_0$mean
p20_2$method <- rep("FullData", nrow(p20_2))
dim(p20m <- rbind(p20_0, p20_1, p20_2))


###---Combine all row-wise
dim(fdat <- as.data.frame(rbind(p100m, p80m, p60m, p40m, p20m)))
fdat$pop_cover <- factor(fdat$pop_cover)
fdat$bld_cover <- factor(fdat$bld_cover)
fdat$method <- factor(fdat$method)

table(fdat$pop_cover, fdat$method)

##--save
write.csv(p100m, "cu_survey_data_100_percent_coverage.csv")
write.csv(p80m, "cu_survey_data_80_percent_coverage.csv")
write.csv(p60m, "cu_survey_data_60_percent_coverage.csv")
write.csv(p40m, "cu_survey_data_40_percent_coverage.csv")
write.csv(p20m, "cu_survey_data_20_percent_coverage.csv")



table(p100m$method)

require(ggpubr)


#####################
###
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

fdat$pop_cover2 <- factor(fdat$pop_cover, levels=var_names)


variable_names <- list(
  "100" = "100% \n Survey Coverage",
  "80" = "80% \n Survey Coverage",
  "60" = "60% \n Survey Coverage",
  "40" = "40% \n Survey Coverage",
  "20" = "20% \n Survey Coverage" 
)

grp <- levels(fdat$method)

variable_labeller <- function(variable,value){
  if (variable=='pop_cover2') {
    return(variable_names[value])
  } else {
    return(grp)
  }
}


#--------------------------plot all survey coverage and all satellite coverage props together


ps <- fdat %>%
  ggplot(aes(x=pop, y=mean))+
    geom_point(data=.%>% filter(bld_cover==100))+
  geom_smooth(aes(colour=bld_cover),method="lm", se=F)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        #panel.border = element_blank(), 
        panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_grid(pop_cover2~method, scales="free", 
             space="free_x",  labeller= variable_labeller)
ps

plots <- ggpar(ps, xlab="Observed Population(Count)", ylab="Predicted Population(Count)",
               legend = "right", legend.title = "Satellite Observation \n Coverage (%)",
               palette = c("#00AFBB", "magenta", "#FC4E07", "#0D0887FF", "#993333", "#2CA25F","#FFE200","red"),
               #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
               #[8] "#FFC600" "#FFE200" "#FFFF00"
               #"#E5F5F9" "#99D8C9" "#FFE200"
               #colour = "bld_cover",
               #shape= "bld_cover",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(16),
               font.y = c(16),
               ylim= c(0, max(p100m$pop)),
               xtickslab.rt = 45, ytickslab.rt = 45)

plots
##---------------------------------------------------------------------------------
##--recode survey coverage levels
levels(fdat$pop_cover2) <- c("Survey \n Coverage:100%", 
                             "Survey \n Coverage:80%",
                             "Survey \n Coverage:60%", 
                             "Survey \n Coverage:40%", 
                             "Survey \n Coverage:20%")

##--Recode Satellite Observations levels
var_namesb <- c(
  "100",
  "95",
  "90",
  "85" ,
  "80",
  "75",
  "70" ,
  "65"
)
fdat$bld_cover2 <- factor(fdat$bld_cover, levels=var_namesb)
table(fdat$bld_cover2)

levels(fdat$bld_cover2)  <- c(
  "Satellite \n Coverage:100%",
  "Satellite \n Coverage:95%",
  "Satellite \n Coverage:90%",
  "Satellite \n Coverage:85%",
  "Satellite \n Coverage:80%", 
  "Satellite \n Coverage:75%",
  "Satellite \n Coverage:70%",
  "Satellite \n Coverage:65%"
)
#===============================================================
p2 <- ggplot(fdat, aes(x=pop, y=mean))+
  geom_point(data=.%>% filter(method=="FullData"))+
  #geom_point()+
  geom_smooth(aes(colour=method),method="lm", se=F)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 11),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        #panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_grid(bld_cover2 ~ pop_cover2, scales="free", 
             space="free_x")
p2

plot2 <- ggpar(p2, xlab="Observed Population(Count)", ylab="Predicted Population(Count)",
               legend = "right", legend.title = "Modelling Method",
               palette = c("#00AFBB","blue", "red"),
               #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
               #[8] "#FFC600" "#FFE200" "#FFFF00"
               #"#E5F5F9" "#99D8C9" "#FFE200"
               #colour = "bld_cover",
               #shape= "bld_cover",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(16),
               font.y = c(16),
               ylim= c(0, max(p100m$pop)),
               xtickslab.rt = 45, ytickslab.rt = 45)

plot2

#=====================================================================================


p3 <- ggplot(fdat, aes(x=pop, y=mean))+
  geom_point(data=.%>% filter(method=="FullData"))+
  #geom_point()+
  geom_smooth(aes(colour=method),method="lm", se=T)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        #panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_grid(bld_cover2 ~ pop_cover2, scales="free", 
             space="free_x")
p3

plot3a <- ggpar(p3, xlab="Observed Population(Count)", ylab="Predicted Population(Count)",
               legend = "right", legend.title = "Modelling Method",
               #palette = c("#00AFBB","blue", "red"),
               #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
               #[8] "#FFC600" "#FFE200" "#FFFF00"
               #"#E5F5F9" "#99D8C9" "#FFE200"
               #colour = "bld_cover",
               #shape= "bld_cover",
               font.label = list(size = 15, face = "bold", color ="red"),
               font.x = c(16),
               font.y = c(16),
               #ylim= c(0, max(p100m$pop)),
               xtickslab.rt = 45, ytickslab.rt = 45)

plot3a
ggsave(plot3a, file="all_data_scatter_3a.tiff", scale=1)
#======

plot3b <- ggpar(p3, xlab="Observed Population(Count)", ylab="Predicted Population(Count)",
                legend = "right", legend.title = "Modelling Method",
                #palette = c("#00AFBB","blue", "red"),
                #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
                #[8] "#FFC600" "#FFE200" "#FFFF00"
                #"#E5F5F9" "#99D8C9" "#FFE200"
                #colour = "bld_cover",
                #shape= "bld_cover",
                font.label = list(size = 13, face = "bold", color ="red"),
                font.x = c(14),
                font.y = c(14),
                ylim= c(0, max(p100m$pop)),
                xtickslab.rt = 45, ytickslab.rt = 45)

plot3b


#ggsave("mtcars.pdf", width = 20, height = 20, units = "cm")

ggsave(plot3b, file="all_data_scatter_3b.tiff", scale=1)

#============================

###-----Boxplots

pvio <- ggplot(fdat, aes(bld_cover, mean, colour=method))+ 
  geom_boxplot(notch=TRUE)+
  #geom_hline(yintercept=mean(fdat$pop[fdat$method=="FullData"]), linetype="dashed", color = "red")+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_grid(pop_cover2~method, scales="free")


plotvio <- ggpar(pvio, xlab="Satellite Coverage(%)", ylab="Predicted Population(Count)",
                legend = "right", legend.title = "Modelling Method",
                #palette = c("#00AFBB","blue", "red"),
                #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
                #[8] "#FFC600" "#FFE200" "#FFFF00"
                #"#E5F5F9" "#99D8C9" "#FFE200"
                #colour = "bld_cover",
                #shape= "bld_cover",
                font.label = list(size = 15, face = "bold", color ="red"),
                font.x = c(16),
                font.y = c(16),
                ylim= c(0, max(p100m$pop)),
                xtickslab.rt = 45, ytickslab.rt = 45)

plotvio
plotvio2 <- plotvio + 
  geom_hline(yintercept=median(fdat$pop[fdat$method=="FullData"]), linetype="dashed", color = "black", size=1)
plotvio2
ggsave(plotvio2, file="all_data_boxplot.tiff", scale=1)


###-----------------------------------------------------------------------------------
####============------------------------
#---------------National estimates------------
dat_nat0<- list()
dat_nat1 <- list()
for(j in 1:n.pop.cover)
{
  pathp <- paste0(out_path,"/outputs_for_", pop.cover[j],"%","_pop_count")
  for(k in 1:n.bld.cover)
  {
    pathb <- paste0(pathp,"/", bld.cover[k], "%","_bldg_count")
    nat0 <- read.csv(paste0(pathb, "/national_estimates_0.csv"))
    nat1 <- read.csv(paste0(pathb, "/national_estimates_1.csv"))
    names(nat1); names(nat0)
    nat0$estimates1 <- nat1$estimates
    
    nat0$bld_cover <- rep(bld.cover[k], nrow(nat0))
    nat0$bld_cover1 <- nat0$bld_cover
    nat0$pop_cover <- rep(pop.cover[j], nrow(nat0))
    nat0$pop_cover1 <- nat0$pop_cover
    #nat0$method <- rep(c("onestep", "twosteps"),nrow(nat0)/2)
    dat_nat0[[k]] = nat0
  }
  dat_nat1[[j]] = dat_nat0
} 
#dat_nat1  #- very huge and 

##--Extract 
ndatpop100 <- dat_nat1[[1]]
ndatpop80 <- dat_nat1[[2]]
ndatpop60 <- dat_nat1[[3]]
ndatpop40 <- dat_nat1[[4]]
ndatpop20 <- dat_nat1[[5]]

class(ndatpop100)
names(ndatpop100)
##---Select variables to include
Var2Include2 <- c("estimates", "estimates1", "bld_cover", "bld_cover1", "pop_cover", "pop_cover1")


##########
##----Extract data for 100% survey coverage
dim(np100 <- rbind(ndatpop100[[1]][,Var2Include2],
                   ndatpop100[[2]][,Var2Include2],
                   ndatpop100[[3]][,Var2Include2],
                   ndatpop100[[4]][,Var2Include2],
                   ndatpop100[[5]][,Var2Include2],
                   ndatpop100[[6]][,Var2Include2],
                   ndatpop100[[7]][,Var2Include2],
                   ndatpop100[[8]][,Var2Include2]))

np100$measure <- rep(c("mean", "lower", "median", "upper"), nrow(np100)/4)


names(np100)
dim(np100_0 <- np100[,names(np100)[c(1,3,5,7)]])
np100_0$method <- rep("BHM", nrow(np100_0))
np100_0$method  <- factor(np100_0$method)
np100_0$bld_cover  <- factor(np100_0$bld_cover)

####--
dim(np100_1 <- np100[,names(np100)[c(2,4,6,7)]])
np100_1$method <- rep("TSBHM", nrow(np100_1))
np100_1$method  <- factor(np100_1$method)
np100_1$bld_cover1  <- factor(np100_1$bld_cover1)

names(np100_1)= names(np100_0)

##----Extract data for 80% survey coverage
dim(np80 <- rbind(ndatpop80[[1]][,Var2Include2],
                  ndatpop80[[2]][,Var2Include2],
                  ndatpop80[[3]][,Var2Include2],
                  ndatpop80[[4]][,Var2Include2],
                  ndatpop80[[5]][,Var2Include2],
                  ndatpop80[[6]][,Var2Include2],
                  ndatpop80[[7]][,Var2Include2],
                  ndatpop80[[8]][,Var2Include2]))

np80$measure <- rep(c("mean", "lower", "median", "upper"), nrow(np80)/4)


names(np80)
dim(np80_0 <- np80[,names(np80)[c(1,3,5,7)]])
np80_0$method <- rep("BHM", nrow(np80_0))
np80_0$method  <- factor(np80_0$method)
np80_0$bld_cover  <- factor(np80_0$bld_cover)

####--extract for TSBHM
dim(np80_1 <- np80[,names(np80)[c(2,4,6,7)]])
np80_1$method <- rep("TSBHM", nrow(np80_1))
np80_1$method  <- factor(np80_1$method)
np80_1$bld_cover1  <- factor(np80_1$bld_cover1)

names(np80_1)= names(np80_0)



##----Extract data for 60% survey coverage
dim(np60 <- rbind(ndatpop60[[1]][,Var2Include2],
                  ndatpop60[[2]][,Var2Include2],
                  ndatpop60[[3]][,Var2Include2],
                  ndatpop60[[4]][,Var2Include2],
                  ndatpop60[[5]][,Var2Include2],
                  ndatpop60[[6]][,Var2Include2],
                  ndatpop60[[7]][,Var2Include2],
                  ndatpop60[[8]][,Var2Include2]))

np60$measure <- rep(c("mean", "lower", "median", "upper"), nrow(np60)/4)


names(np60)
dim(np60_0 <- np60[,names(np60)[c(1,3,5,7)]])
np60_0$method <- rep("BHM", nrow(np60_0))
np60_0$method  <- factor(np60_0$method)
np60_0$bld_cover  <- factor(np60_0$bld_cover)

####--
dim(np60_1 <- np60[,names(np60)[c(2,4,6,7)]])
np60_1$method <- rep("TSBHM", nrow(np60_1))
np60_1$method  <- factor(np60_1$method)
np60_1$bld_cover1  <- factor(np60_1$bld_cover1)

names(np60_1)= names(np60_0)


##----Extract data for 40% survey coverage
dim(np40 <- rbind(ndatpop40[[1]][,Var2Include2],
                  ndatpop40[[2]][,Var2Include2],
                  ndatpop40[[3]][,Var2Include2],
                  ndatpop40[[4]][,Var2Include2],
                  ndatpop40[[5]][,Var2Include2],
                  ndatpop40[[6]][,Var2Include2],
                  ndatpop40[[7]][,Var2Include2],
                  ndatpop40[[8]][,Var2Include2]))

np40$measure <- rep(c("mean", "lower", "median", "upper"), nrow(np40)/4)


names(np40)
dim(np40_0 <- np40[,names(np40)[c(1,3,5,7)]])
np40_0$method <- rep("BHM", nrow(np40_0))
np40_0$method  <- factor(np40_0$method)
np40_0$bld_cover  <- factor(np40_0$bld_cover)

####--
dim(np40_1 <- np40[,names(np40)[c(2,4,6,7)]])
np40_1$method <- rep("TSBHM", nrow(np40_1))
np40_1$method  <- factor(np40_1$method)
np40_1$bld_cover1  <- factor(np40_1$bld_cover1)

names(np40_1)= names(np40_0)


##----Extract data for 20% survey coverage
dim(np20 <- rbind(ndatpop20[[1]][,Var2Include2],
                  ndatpop20[[2]][,Var2Include2],
                  ndatpop20[[3]][,Var2Include2],
                  ndatpop20[[4]][,Var2Include2],
                  ndatpop20[[5]][,Var2Include2],
                  ndatpop20[[6]][,Var2Include2],
                  ndatpop20[[7]][,Var2Include2],
                  ndatpop20[[8]][,Var2Include2]))

np20$measure <- rep(c("mean", "lower", "median", "upper"), nrow(np20)/4)


names(np20)
dim(np20_0 <- np20[,names(np20)[c(1,3,5,7)]])
np20_0$method <- rep("BHM", nrow(np20_0))
np20_0$method  <- factor(np20_0$method)
np20_0$bld_cover  <- factor(np20_0$bld_cover)

####--
dim(np20_1 <- np20[,names(np20)[c(2,4,6,7)]])
np20_1$method <- rep("TSBHM", nrow(np20_1))
np20_1$method  <- factor(np20_1$method)
np20_1$bld_cover1  <- factor(np20_1$bld_cover1)

names(np20_1)= names(np20_0)


###---join all ndata
dim(np100m <- rbind(np100_0, np100_1))
dim(np80m <- rbind(np80_0, np80_1))
dim(np60m <- rbind(np60_0, np60_1))
dim(np40m <- rbind(np40_0, np40_1))
dim(np20m <- rbind(np20_0, np20_1))

##----
dim(nfdat <- as.data.frame(rbind(np100m, np80m, np60m, np40m, np20m)))
nfdat$pop_cover <- factor(nfdat$pop_cover)
nfdat$bld_cover <- factor(nfdat$bld_cover)
nfdat$method <- factor(nfdat$method)

table(nfdat$pop_cover, nfdat$method)

dim(nfdat_mean <- nfdat[nfdat$measure=="mean",]) #---select only mean
dim(nfdat_lower <- nfdat[nfdat$measure=="lower",]) #---select only lower bound 2.5%
dim(nfdat_median <- nfdat[nfdat$measure=="median",]) #---select only median
dim(nfdat_upper <- nfdat[nfdat$measure=="upper",]) #---select only upper bound 97.5%

##----
names(nfdat_mean) <- c("mean", "bld_cover", "pop_cover", "measure", "method")
nfdat_mean$lower <- nfdat_lower$estimates
nfdat_mean$upper <- nfdat_upper$estimates


###---plots-----------------------
var_names <- c(
  "100",
  "80",
  "60",
  "40" ,
  "20"
)

nfdat_mean$pop_cover2 <- factor(nfdat_mean$pop_cover, levels=var_names)


variable_names <- list(
  "100" = "100% \n Survey Coverage",
  "80" = "80% \n Survey Coverage",
  "60" = "60% \n Survey Coverage",
  "40" = "40% \n Survey Coverage",
  "20" = "20% \n Survey Coverage" 
)

ngrp <- levels(nfdat_mean$method)

variable_labeller <- function(variable,value){
  if (variable=='pop_cover2') {
    return(variable_names[value])
  } else {
    return(ngrp)
  }
}


levels(nfdat_mean$pop_cover2) <- c(
   "Survey Coverage:\n 100%",
  "Survey Coverage: \n 80%",
  "Survey Coverage: \n 60%",
  "Survey Coverage: \n 40%",
  "Survey Coverage: \n 20%" 
)

# Create a simple example dataset

pp <- ggplot(nfdat_mean, aes(bld_cover, mean, colour=method))+ 
  geom_linerange(aes(ymin = lower, ymax = upper), size=1) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size=0.5)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=1,width = 0.5)+
  geom_line(aes(group = method), size=1) +
  geom_hline(yintercept=11625153 , linetype="dashed", color = "black")+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))+
  facet_wrap(~pop_cover2, scales="free", nrow=3)

plotp<- ggpar(pp, xlab="Satellite Coverage(%)", ylab="Predicted Population(Count)",
                legend = "right", legend.title = "Modelling Method",
                #palette = c("#00AFBB","blue", "red"),
                #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
                #[8] "#FFC600" "#FFE200" "#FFFF00"
                #"#E5F5F9" "#99D8C9" "#FFE200"
                #colour = "bld_cover",
                #shape= "bld_cover",
                font.label = list(size = 15, face = "bold", color ="red"),
                font.x = c(16),
                font.y = c(16),
                #ylim= c(0, max(p100m$pop)),
                xtickslab.rt = 45, ytickslab.rt = 45)

plotp


#ggsave("mtcars.pdf", width = 20, height = 20, units = "cm")

ggsave(plotp, file="national_all_data_error_bar.tiff", scale=1)

write.csv(nfdat_mean, "national_all_estimates2.csv")

nfdat_mean$per_error <- ((nfdat_mean$mean/11643074)-1)*100
nfdat_mean$per_error2 <- (nfdat_mean$mean/11643074)*100

nfdat_mean$bias <- ((nfdat_mean$mean-11643074)/11643074)*100
ggplot(nfdat_mean, aes(fill=method, y=per_error2, x=pop_cover)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~bld_cover)


ggplot(nfdat_mean, aes(fill=method, y=bias, x=pop_cover)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~bld_cover)


ggplot(nfdat_mean, aes(color=method, y=bias, x=pop_cover)) + 
  geom_point()+
  geom_line(aes(group = method), size=1) +
  facet_wrap(~bld_cover)


###
ggplot(nfdat_mean, aes(color=bld_cover, y=bias, x=pop_cover)) + 
  geom_point()+
  geom_line(aes(group = bld_cover), size=1) +
  facet_wrap(~method)


##--bias
ggplot(nfdat_mean, aes(color=pop_cover, y=bias, x=bld_cover)) + 
  geom_point()+
  geom_line(aes(group = pop_cover), size=1) +
  facet_wrap(~method)

##
ggplot(nfdat_mean, aes(fill=pop_cover, y=bias, x=bld_cover)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~method)


##
ggplot(nfdat_mean, aes(fill=method, y=bias, x=pop_cover)) + 
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~bld_cover)




###-------
dat_m1 <- nfdat_mean[nfdat_mean$method=="BHM",]
dat_m2 <- nfdat_mean[nfdat_mean$method=="TSBHM",]
names(dat_m2)[11] <- "bias2"
names(dat_m2)[5] <- "method2"

names(nfdat_mean)
dat_m <- as.data.frame(cbind(dat_m1, dat_m2[,c("method2","bias2")]))
names(dat_m)


dat_mm <- dat_m
dat_mm$rbias <- (1-(dat_mm$bias2/dat_mm$bias))*100

dat_mm <- dat_m[dat_m$pop_cover!=100 & dat_m$bld_cover!=100,]
dat_mm$rbias <- (1-(dat_mm$bias2/dat_mm$bias))*100

c(min(dat_mm$rbias), max(dat_mm$rbias))

ggplot(dat_mm, aes( y=rbias, x=bld_cover, col="pop_cover,")) + 
  geom_point()+
  geom_()+
  facet_wrap(~pop_cover)




library(ggplot2)
ggp <- ggplot(perf)  + 
  geom_bar(aes(x=year, y=course),stat="identity", fill="cyan",colour="#006000")+
  geom_line(aes(x=year, y=100*penroll),stat="identity",color="red",size=2)+
  labs(title= "Courses vs Students Enrolled in GeeksforGeeks",
       x="Year",y="Number of Courses Sold")+
  scale_y_continuous(sec.axis=sec_axis(
    ~.*0.01,name="Percentage of Students Enrolled", labels=scales::percent))
ggp



ggplot(dat_mm, aes(color=pop_cover, y=rbias, x=bld_cover)) + 
  geom_point()+
  geom_line( size=1) +
  facet_wrap(~pop_cover)



###--Plot histogram
gghistogram(wdata, x = "weight",
            add = "mean", rug = TRUE,
            color = "sex", fill = "sex",
            palette = c("#00AFBB", "#E7B800"))



density.p <- ggdensity(iris, x = "Sepal.Length", 
                       fill = "Species", palette = "jco")




library(dplyr)
library(tidyr)
pop_1950 <- pop %>%
  gather(Sex, Number, females, males)


ggplot(pop) +
  geom_col(aes(x=age_groups, y=females, fill = "red")) +
  geom_col(aes(x=age_groups, y=-males, fill = "blue")) + 
  scale_y_continuous(breaks = seq(-200000, 200000, 50000), 
                     labels = paste0(as.character(c(seq(200, 0, -50), seq(50, 200, 50))))) + 
  coord_flip() + 
  scale_fill_discrete(name = "Gender", labels = c("Female", "Male"))+
  xlab("Age")+
  ylab("Population (000's)")+
  theme_minimal()




# Dot plot (dp)
dp <- ggdotplot(ToothGrowth, x = "dose", y = "len",
                color = "dose", palette = "jco", binwidth = 1)
dp




ggviolin(ToothGrowth, x = "dose", y = "len", 
         fill = "dose",
         palette = c("#00AFBB", 
                     "#E7B800", "#FC4E07"),
         add = "boxplot", 
         add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, 
                     label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)# Add global the p-value 



dat_mm

dat_mm$names <- factor(paste0(dat_mm$bld_cover, "_for_", dat_mm$pop_cover))

#'#999999','#E69F00','#56B4E9'

##--relative bias
#####--Using numbers
lpp <- ggdotchart(dat_mm, x = "names", y = "rbias",
           color = "pop_cover",                                # Color by groups
           palette = c("#00AFBB", "#E7B800", "#FC4E07", "#999999"),            # jco journal color palett. see ?ggpar
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           size = 1,
           rotate = TRUE,                                # Rotate vertically
           group = "pop_cover",
           #egend.text=element_text(size=22),# Order by groups
           dot.size = 18,                                 # Large dot size
           label = round(dat_mm$rbias),  
           #scale_x_continuous(limits = c(0, 100), breaks = c(20, 40, 60, 80,100)),# Add mpg values as dot labels
           font.label = list(color = "white", size = 22,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)

plot_lpp<- ggpar(lpp, ylab="Perecentage reduction in relative bias", xlab="Sattelite-Survey Percentage Coverages",
              legend = "right", legend.title=element_text("Survey Coverage (%)",size=22) ,
              font.legend=c(18),
              
             # legend.text=element_text(size=22),
              #main ="Lollipop plot of percentage relative  bias reduction by TSBHM",
              #palette = c("#00AFBB","blue", "red"),
              #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
              #[8] "#FFC600" "#FFE200" "#FFFF00"
              #"#E5F5F9" "#99D8C9" "#FFE200"
              #colour = "bld_cover",
              #shape= "bld_cover",
             #scale_x_continuous(limits=c(0,100), breaks=c(20, 40,60, 80, 100)),
              font.label = list(size = 18, face = "bold", color ="red"),
              font.x = c(22),
              font.y = c(18),
              font.main=c(14),
              font.xtickslab =c(18),
              font.ytickslab =c(16),
              #xticks.by = TRUE,
              #ylim= c(0, max(p100m$pop)),
              xtickslab.rt = 45, ytickslab.rt = 45)

plot_lpp



#############
###----OBTAIN OBSERVED TOTAL COUNT FOR EACH SURVEY COVERAGE
p100_summary <- p100_0 %>% group_by(bld_cover) %>%
  summarise(popt=sum(popm, na.rm=T))
p100_summary 

p80_summary <- p80_0 %>% group_by(bld_cover) %>%
  summarise(popt=sum(popm, na.rm=T))
p80_summary 

p60_summary <- p60_0 %>% group_by(bld_cover) %>%
  summarise(popt=sum(popm, na.rm=T))
p60_summary 

p40_summary <- p40_0 %>% group_by(bld_cover) %>%
  summarise(popt=sum(popm, na.rm=T))
p40_summary 

p20_summary <- p20_0 %>% group_by(bld_cover) %>%
  summarise(popt=sum(popm, na.rm=T))
p20_summary 

##--save

#################-----MAPS---------------------------------------------------------------
###-----------------------------------------------------------------------------
#####---------------------------------------------------------------------------------
# for loading our data
library(raster);library(readr);library(readxl);library(sf)
# for datasets
library(maps);library(spData)
# for creating animations
library(magick)
# for plotting
library(grid);library(tmap);library(viridis); library(tmaptools)

##=================================================================================================
###-----Simulate the coordinates based on the lon and lat of PNG
shp <- readOGR(dsn = paste0(path), layer = "PNG_CU_32100_B") #combined shapefile 
plot(shp)



#post_maps <- function(cu.est, shp)
#{

shp.cu <- shp
shp.cu$obs1 <- p100_0_b1$popm
shp.cu$mean10 <- p100_0_b1$mean
shp.cu$mean11<- p100_1_b1$mean
shp.cu$uncert10 <- (p100_0_b1$upper-p100_0_b1$lower)/p100_0_b1$mean
shp.cu$uncert11 <- (p100_1_b1$upper-p100_1_b1$lower)/p100_1_b1$mean


shp.cu$obs2 <- p80_0_b1$popm
shp.cu$mean20 <- p80_0_b1$mean
shp.cu$mean21 <- p80_1_b1$mean
shp.cu$uncert20 <- (p80_0_b1$upper-p80_0_b1$lower)/p80_0_b1$mean
shp.cu$uncert21 <- (p80_1_b1$upper-p80_1_b1$lower)/p80_1_b1$mean


shp.cu$obs3 <- p60_0_b1$popm
shp.cu$mean30 <- p60_0_b1$mean
shp.cu$mean31 <- p60_1_b1$mean
shp.cu$uncert30 <- (p60_0_b1$upper-p60_0_b1$lower)/p60_0_b1$mean
shp.cu$uncert31 <- (p60_1_b1$upper-p60_1_b1$lower)/p60_1_b1$mean


shp.cu$obs4 <- p40_0_b1$popm
shp.cu$mean40 <- p40_0_b1$mean
shp.cu$mean41 <- p40_1_b1$mean
shp.cu$uncert40 <- (p40_0_b1$upper-p40_0_b1$lower)/p40_0_b1$mean
shp.cu$uncert41 <- (p40_1_b1$upper-p40_1_b1$lower)/p40_1_b1$mean


shp.cu$obs5 <- p20_0_b1$popm
shp.cu$mean50 <- p20_0_b1$mean
shp.cu$mean51 <- p20_1_b1$mean
shp.cu$uncert50 <- (p20_0_b1$upper-p20_0_b1$lower)/p20_0_b1$mean
shp.cu$uncert51 <- (p20_1_b1$upper-p20_1_b1$lower)/p20_1_b1$mean



###
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM = spTransform(shp.cu, crs.UTM)

names(shpcu_UTM)
#plot(shpcu_UTM)
#spplot(shpcu_UTM, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
#max(shp100 <- c(shpcu_UTM$obs1, shpcu_UTM$obs2, shpcu_UTM$obs3, shpcu_UTM$obs4, shpcu_UTM$obs5,
 #               shpcu_UTM$pred10, shpcu_UTM$pred20, shpcu_UTM$pred30, shpcu_UTM$pred40, shpcu_UTM$pred50,
#                shpcu_UTM$pred11, shpcu_UTM$pred21, shpcu_UTM$pred31, shpcu_UTM$pred41, shpcu_UTM$pred51), na.rm=T)


#max(shp100u <- c(shpcu_UTM$uncert10, shpcu_UTM$uncert20, shpcu_UTM$uncert30, shpcu_UTM$uncert40, shpcu_UTM$uncert50,
 #               shpcu_UTM$uncert11, shpcu_UTM$uncert21, shpcu_UTM$uncert31, shpcu_UTM$uncert41, shpcu_UTM$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
 # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.005,0.01,0.015,0.030, 0.050,0.075))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



###########################################################################
#####--------------95% Satellite Observation coverage 
############################################################################
shp.cu2 <- shp
shp.cu2$obs1 <- p100_0_b2$popm
shp.cu2$mean10 <- p100_0_b2$mean
shp.cu2$mean11<- p100_1_b2$mean
shp.cu2$uncert10 <- (p100_0_b2$upper-p100_0_b2$lower)/p100_0_b2$mean
shp.cu2$uncert11 <- (p100_1_b2$upper-p100_1_b2$lower)/p100_1_b2$mean


shp.cu2$obs2 <- p80_0_b2$popm
shp.cu2$mean20 <- p80_0_b2$mean
shp.cu2$mean21 <- p80_1_b2$mean
shp.cu2$uncert20 <- (p80_0_b2$upper-p80_0_b2$lower)/p80_0_b2$mean
shp.cu2$uncert21 <- (p80_1_b2$upper-p80_1_b2$lower)/p80_1_b2$mean


shp.cu2$obs3 <- p60_0_b2$popm
shp.cu2$mean30 <- p60_0_b2$mean
shp.cu2$mean31 <- p60_1_b2$mean
shp.cu2$uncert30 <- (p60_0_b2$upper-p60_0_b2$lower)/p60_0_b2$mean
shp.cu2$uncert31 <- (p60_1_b2$upper-p60_1_b2$lower)/p60_1_b2$mean


shp.cu2$obs4 <- p40_0_b2$popm
shp.cu2$mean40 <- p40_0_b2$mean
shp.cu2$mean41 <- p40_1_b2$mean
shp.cu2$uncert40 <- (p40_0_b2$upper-p40_0_b2$lower)/p40_0_b2$mean
shp.cu2$uncert41 <- (p40_1_b2$upper-p40_1_b2$lower)/p40_1_b2$mean


shp.cu2$obs5 <- p20_0_b2$popm
shp.cu2$mean50 <- p20_0_b2$mean
shp.cu2$mean51 <- p20_1_b2$mean
shp.cu2$uncert50 <- (p20_0_b2$upper-p20_0_b2$lower)/p20_0_b2$mean
shp.cu2$uncert51 <- (p20_1_b2$upper-p20_1_b2$lower)/p20_1_b2$mean



###
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM2 = spTransform(shp.cu2, crs.UTM)

names(shpcu_UTM2)
#plot(shpcu_UTM2)
#spplot(shpcu_UTM2, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM2$obs1, shpcu_UTM2$obs2, shpcu_UTM2$obs3, shpcu_UTM2$obs4, shpcu_UTM2$obs5,
                shpcu_UTM2$pred10, shpcu_UTM2$pred20, shpcu_UTM2$pred30, shpcu_UTM2$pred40, shpcu_UTM2$pred50,
                shpcu_UTM2$pred11, shpcu_UTM2$pred21, shpcu_UTM2$pred31, shpcu_UTM2$pred41, shpcu_UTM2$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM2$uncert10, shpcu_UTM2$uncert20, shpcu_UTM2$uncert30, shpcu_UTM2$uncert40, shpcu_UTM2$uncert50,
                 shpcu_UTM2$uncert11, shpcu_UTM2$uncert21, shpcu_UTM2$uncert31, shpcu_UTM2$uncert41, shpcu_UTM2$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM2)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.3,0.4, 0.6,0.95))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


###########################################################################
#####--------------90% Satellite Observation coverage 
############################################################################
shp.cu3 <- shp
shp.cu3$obs1 <- p100_0_b3$popm
shp.cu3$mean10 <- p100_0_b3$mean
shp.cu3$mean11<- p100_1_b3$mean
shp.cu3$uncert10 <- (p100_0_b3$upper-p100_0_b3$lower)/p100_0_b3$mean
shp.cu3$uncert11 <- (p100_1_b3$upper-p100_1_b3$lower)/p100_1_b3$mean


shp.cu3$obs2 <- p80_0_b3$popm
shp.cu3$mean20 <- p80_0_b3$mean
shp.cu3$mean21 <- p80_1_b3$mean
shp.cu3$uncert20 <- (p80_0_b3$upper-p80_0_b3$lower)/p80_0_b3$mean
shp.cu3$uncert21 <- (p80_1_b3$upper-p80_1_b3$lower)/p80_1_b3$mean


shp.cu3$obs3 <- p60_0_b3$popm
shp.cu3$mean30 <- p60_0_b3$mean
shp.cu3$mean31 <- p60_1_b3$mean
shp.cu3$uncert30 <- (p60_0_b3$upper-p60_0_b3$lower)/p60_0_b3$mean
shp.cu3$uncert31 <- (p60_1_b3$upper-p60_1_b3$lower)/p60_1_b3$mean


shp.cu3$obs4 <- p40_0_b3$popm
shp.cu3$mean40 <- p40_0_b3$mean
shp.cu3$mean41 <- p40_1_b3$mean
shp.cu3$uncert40 <- (p40_0_b3$upper-p40_0_b3$lower)/p40_0_b3$mean
shp.cu3$uncert41 <- (p40_1_b3$upper-p40_1_b3$lower)/p40_1_b3$mean


shp.cu3$obs5 <- p20_0_b3$popm
shp.cu3$mean50 <- p20_0_b3$mean
shp.cu3$mean51 <- p20_1_b3$mean
shp.cu3$uncert50 <- (p20_0_b3$upper-p20_0_b3$lower)/p20_0_b3$mean
shp.cu3$uncert51 <- (p20_1_b3$upper-p20_1_b3$lower)/p20_1_b3$mean



###
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM3 = spTransform(shp.cu3, crs.UTM)

names(shpcu_UTM3)
#plot(shpcu_UTM3)
#spplot(shpcu_UTM3, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM3$obs1, shpcu_UTM3$obs2, shpcu_UTM3$obs3, shpcu_UTM3$obs4, shpcu_UTM3$obs5,
                shpcu_UTM3$pred10, shpcu_UTM3$pred20, shpcu_UTM3$pred30, shpcu_UTM3$pred40, shpcu_UTM3$pred50,
                shpcu_UTM3$pred11, shpcu_UTM3$pred21, shpcu_UTM3$pred31, shpcu_UTM3$pred41, shpcu_UTM3$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM3$uncert10, shpcu_UTM3$uncert20, shpcu_UTM3$uncert30, shpcu_UTM3$uncert40, shpcu_UTM3$uncert50,
                 shpcu_UTM3$uncert11, shpcu_UTM3$uncert21, shpcu_UTM3$uncert31, shpcu_UTM3$uncert41, shpcu_UTM3$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20)
  )+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM3)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 




###########################################################################
#####--------------85% Satellite Observation coverage 
############################################################################
shp.cu4 <- shp
shp.cu4$obs1 <- p100_0_b4$popm
shp.cu4$mean10 <- p100_0_b4$mean
shp.cu4$mean11<- p100_1_b4$mean
shp.cu4$uncert10 <- (p100_0_b4$upper-p100_0_b4$lower)/p100_0_b4$mean
shp.cu4$uncert11 <- (p100_1_b4$upper-p100_1_b4$lower)/p100_1_b4$mean


shp.cu4$obs2 <- p80_0_b4$popm
shp.cu4$mean20 <- p80_0_b4$mean
shp.cu4$mean21 <- p80_1_b4$mean
shp.cu4$uncert20 <- (p80_0_b4$upper-p80_0_b4$lower)/p80_0_b4$mean
shp.cu4$uncert21 <- (p80_1_b4$upper-p80_1_b4$lower)/p80_1_b4$mean


shp.cu4$obs3 <- p60_0_b4$popm
shp.cu4$mean30 <- p60_0_b4$mean
shp.cu4$mean31 <- p60_1_b4$mean
shp.cu4$uncert30 <- (p60_0_b4$upper-p60_0_b4$lower)/p60_0_b4$mean
shp.cu4$uncert31 <- (p60_1_b4$upper-p60_1_b4$lower)/p60_1_b4$mean


shp.cu4$obs4 <- p40_0_b4$popm
shp.cu4$mean40 <- p40_0_b4$mean
shp.cu4$mean41 <- p40_1_b4$mean
shp.cu4$uncert40 <- (p40_0_b4$upper-p40_0_b4$lower)/p40_0_b4$mean
shp.cu4$uncert41 <- (p40_1_b4$upper-p40_1_b4$lower)/p40_1_b4$mean


shp.cu4$obs5 <- p20_0_b4$popm
shp.cu4$mean50 <- p20_0_b4$mean
shp.cu4$mean51 <- p20_1_b4$mean
shp.cu4$uncert50 <- (p20_0_b4$upper-p20_0_b4$lower)/p20_0_b4$mean
shp.cu4$uncert51 <- (p20_1_b4$upper-p20_1_b4$lower)/p20_1_b4$mean



###
require(rgdal)
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM4 = spTransform(shp.cu4, crs.UTM)

names(shpcu_UTM4)
#plot(shpcu_UTM4)
#spplot(shpcu_UTM4, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM4$obs1, shpcu_UTM4$obs2, shpcu_UTM4$obs3, shpcu_UTM4$obs4, shpcu_UTM4$obs5,
                shpcu_UTM4$pred10, shpcu_UTM4$pred20, shpcu_UTM4$pred30, shpcu_UTM4$pred40, shpcu_UTM4$pred50,
                shpcu_UTM4$pred11, shpcu_UTM4$pred21, shpcu_UTM4$pred31, shpcu_UTM4$pred41, shpcu_UTM4$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM4$uncert10, shpcu_UTM4$uncert20, shpcu_UTM4$uncert30, shpcu_UTM4$uncert40, shpcu_UTM4$uncert50,
                 shpcu_UTM4$uncert11, shpcu_UTM4$uncert21, shpcu_UTM4$uncert31, shpcu_UTM4$uncert41, shpcu_UTM4$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25)
  )+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.25))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM4)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.8,1.20))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 




###########################################################################
#####--------------80% Satellite Observation coverage 
############################################################################
shp.cu5 <- shp
shp.cu5$obs1 <- p100_0_b5$popm
shp.cu5$mean10 <- p100_0_b5$mean
shp.cu5$mean11<- p100_1_b5$mean
shp.cu5$uncert10 <- (p100_0_b5$upper-p100_0_b5$lower)/p100_0_b5$mean
shp.cu5$uncert11 <- (p100_1_b5$upper-p100_1_b5$lower)/p100_1_b5$mean


shp.cu5$obs2 <- p80_0_b5$popm
shp.cu5$mean20 <- p80_0_b5$mean
shp.cu5$mean21 <- p80_1_b5$mean
shp.cu5$uncert20 <- (p80_0_b5$upper-p80_0_b5$lower)/p80_0_b5$mean
shp.cu5$uncert21 <- (p80_1_b5$upper-p80_1_b5$lower)/p80_1_b5$mean


shp.cu5$obs3 <- p60_0_b5$popm
shp.cu5$mean30 <- p60_0_b5$mean
shp.cu5$mean31 <- p60_1_b5$mean
shp.cu5$uncert30 <- (p60_0_b5$upper-p60_0_b5$lower)/p60_0_b5$mean
shp.cu5$uncert31 <- (p60_1_b5$upper-p60_1_b5$lower)/p60_1_b5$mean


shp.cu5$obs4 <- p40_0_b5$popm
shp.cu5$mean40 <- p40_0_b5$mean
shp.cu5$mean41 <- p40_1_b5$mean
shp.cu5$uncert40 <- (p40_0_b5$upper-p40_0_b5$lower)/p40_0_b5$mean
shp.cu5$uncert41 <- (p40_1_b5$upper-p40_1_b5$lower)/p40_1_b5$mean


shp.cu5$obs5 <- p20_0_b5$popm
shp.cu5$mean50 <- p20_0_b5$mean
shp.cu5$mean51 <- p20_1_b5$mean
shp.cu5$uncert50 <- (p20_0_b5$upper-p20_0_b5$lower)/p20_0_b5$mean
shp.cu5$uncert51 <- (p20_1_b5$upper-p20_1_b5$lower)/p20_1_b5$mean



###
require(rgdal)
crs.UTM = CRS("+proj=utm +zone=55 +datum=WGS84 +units=m +no_defs")#png is zone 55
shpcu_UTM5 = spTransform(shp.cu5, crs.UTM)

names(shpcu_UTM5)
#plot(shpcu_UTM5)
#spplot(shpcu_UTM5, "mean")



tmap_options(check.and.fix = TRUE)
tmap_mode("plot")


#
max(shp100 <- c(shpcu_UTM5$obs1, shpcu_UTM5$obs2, shpcu_UTM5$obs3, shpcu_UTM5$obs4, shpcu_UTM5$obs5,
                shpcu_UTM5$pred10, shpcu_UTM5$pred20, shpcu_UTM5$pred30, shpcu_UTM5$pred40, shpcu_UTM5$pred50,
                shpcu_UTM5$pred11, shpcu_UTM5$pred21, shpcu_UTM5$pred31, shpcu_UTM5$pred41, shpcu_UTM5$pred51), na.rm=T)


max(shp100u <- c(shpcu_UTM5$uncert10, shpcu_UTM5$uncert20, shpcu_UTM5$uncert30, shpcu_UTM5$uncert40, shpcu_UTM5$uncert50,
                 shpcu_UTM5$uncert11, shpcu_UTM5$uncert21, shpcu_UTM5$uncert31, shpcu_UTM5$uncert41, shpcu_UTM5$uncert51), na.rm=T)

##=============================================================
#--Observed counts
obs1 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs1", title="Simulated Count",
              #style="fixed", palette="viridis")+
              style="fixed", 
              #as.count=T,
              #n=6,
              #legend.hist=T,
              palette="viridis", 
              legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside =T, legend.text.size=1.5, legend.title.size=1.5)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1.5, size=1, breaks=c(0, 200, 400, 600))+
  #tm_layout(main.title = "Simulated Counts at 100% Survey Coverage") 
  tm_layout(main.title = "") 

##--predicted counts 
pred10 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean10", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = T, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert10 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert10", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4)
  )+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred11 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean11", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert11<-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert11", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------
#--Observed counts
obs2 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs2", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = F,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred20 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean20", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert20 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert20", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred21 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean21", title="Predicted Count",
              style="fixed", legend.hist=F,, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert21 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert21", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--------------------------------------------------------------------------------------------------------------
#--Observed counts
obs3 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs3", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred30 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean30", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert30 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert30", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred31 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean31", title="Predicted Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  # tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert31 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert31", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

#-----------------------------------------------------------------------------------------------------
#--Observed counts
obs4 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs4", title="Observed Count",
              style="fixed", palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred40 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean40", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



##--uncertainty
uncert40 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert40", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


pred41 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean41", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert41<-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert41", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##============================================================================= 
#--Observed counts
obs5 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="obs5", title="Observed Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 

##--predicted counts 
pred50 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean50", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


##--uncertainty
uncert50 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert50", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



pred51 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="mean51", title="Predicted Count",
              style="fixed", legend.hist=F, palette="viridis", legend.show = FALSE,
              breaks=c(1,50,150,300,450, 700,1000,3000))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  #tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 


uncert51 <-  tm_shape(shpcu_UTM5)+ 
  tm_polygons(col="uncert51", title="Uncertainty",
              #style="fixed", 
              legend.hist=F, 
              palette=viridis(n=100, option="C", direction=-1), 
              legend.show =F,
              breaks=c(0,0.1,0.2,0.4,0.6, 0.9,1.4))+
  tm_layout(legend.outside = F, legend.text.size=1, legend.title.size=1)+
  #tm_compass(position = c("right", "bottom"), text.size=1.5)+
  # tm_scale_bar(position = c("left", "bottom"), text.size=1, size=1, breaks=c(0, 200, 400, 600))+
  tm_layout(main.title = "") 



####################################################################################################
######--------------------Perecentage error---------------------------------------------------------
dim(p100mm <- p100m[p100m$method!="FullData",])
p100mm$p_error <- (p100mm$mean/p100mm$pop)*100

dim(p80mm <- p80m[p80m$method!="FullData",])
p80mm$p_error <- (p80mm$mean/p80mm$pop)*100

dim(p60mm <- p60m[p60m$method!="FullData",])
p60mm$p_error <- (p60mm$mean/p60mm$pop)*100

dim(p40mm <- p40m[p40m$method!="FullData",])
p40mm$p_error <- (p40mm$mean/p40mm$pop)*100

dim(p20mm <- p20m[p20m$method!="FullData",])
p20mm$p_error <- (p20mm$mean/p20mm$pop)*100


table(fdat$method)
dim(fdatm <- fdat[fdat$method!="FullData",])
fdatm$R_error <- (fdatm$mean/fdatm$pop)/(sum(fdatm$mean/fdatm$pop))

# plot
bar_plot <- ggplot(fdatm, aes(fill=method, y=R_error, x=pop_cover)) + 
  geom_bar(position="stack", stat="identity")+
  #geom_hline(yintercept=1 , linetype="dashed", color = "black")+
  #scale_fill_viridis(discrete=TRUE, name="") +
  theme_bw()+
  #theme_ipsum() +
  theme(strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11),
        axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14),
        #panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+
  facet_wrap(~bld_cover2, ncol=2)
bar_plot<- ggpar(bar_plot, xlab="Survey Coverage(%)", ylab="Normalized Relative Error Rate",
              legend = "right", legend.title = "Modelling Method",
              #palette = c("#00AFBB","blue", "red"),
              #[1] "#FF0000" "#FF1C00" "#FF3800" "#FF5500" "#FF7100" "#FF8D00" "#FFAA00", "#E5F5F9","#99D8C9" 
              #[8] "#FFC600" "#FFE200" "#FFFF00"
              #"#E5F5F9" "#99D8C9" "#FFE200"
              #colour = "bld_cover",
              #shape= "bld_cover",
              font.label = list(size = 15, face = "bold", color ="red"),
              font.x = c(16),
              font.y = c(16),
              #ylim= c(0, max(p100m$pop)),
              xtickslab.rt = 45, ytickslab.rt = 45)

bar_plot

ggsave(plot3a, file="all_data_scatter_3a.tiff", scale=1)




















