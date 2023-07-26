####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR PAPUA NEW GUNEA BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON 
####--DATE:19 DECEMBER 2022.
###=====================================================================================================

rm(list=ls()) #---Clear workspace
#set.seed(3125)
##set.seed(312)
##---List of potentially key packages
packages <- c("raster", "haven", "sf","sp", "tidyverse","rgdal",
              "lattice", "gridExtra", "devtools", "rlang", "spdep", "viridis", "MASS")

##----Install packages listed and not yet installed
if(length(setdiff(packages, rownames(installed.packages()))) > 0) { 
  install.packages(setdiff(packages, rownames(installed.packages())), type="binary") }

#Installing INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}

#install.packages("INLA")
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries


##-------------------------------------------------------------------------------
###--Specify various file paths (set input and output paths)
#--------------------------------------------------------------------------------
data_path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP000008_UNFPA_PNG/DataOut/ModelOutputs_final/"
results_path <- paste0(data_path, "Results/")

dir.create(paste0(data_path,"Results"), recursive = T, showWarnings = F)
dir.create(paste0(data_path,"Results/plots"), recursive = T, showWarnings = F)
dir.create(paste0(data_path,"Results/updated2"), recursive = T, showWarnings = F)

#------------------------Load datasets-----------------------------------------        
#--------------------------------------------------------------------------------
shp <- readOGR(dsn = paste0(data_path,"gg_model_input_data/Boundary_files/CU"), layer = "PNG_CU_32100_B") #combined shapefile 
#plot(shp);dim(shp); names(shp)

##load Rdata--
load(paste0(data_path,"gg_model_input_data/Survey_data/ModelDataPNG_BI.RData")) #-


###---Explore  and clean data further wherever necessary----------------------------------------
#--------------------------------------------------------------------------
ls()
names(covs)
covDF <- as.data.frame(covs)
#plot(covDF$lon, covDF$lat)
#names(covDF); str(covDF); dim(covDF); dim(shp)# shp has one extra row

shp <- shp[-32100,] #--remove row 32100 from the shapefile to match with the survey data
dim(covDF); dim(shp) #---confirm changes made
#--------------------------------------------------------------------------

##----Extract the coordinates - the centroids of the CUs
#--------------------------------------------------------------------
library(sp)
covDF$lon <- coordinates(shp)[,1]#--add lon-lat to the data
covDF$lat <- coordinates(shp)[,2]
covDF2 <- covDF
#---------------------------------------------------------------------------



#------------------------------------------------------------------------
#---Define some key functions to be used later
#--------------------------------------------------------------------------
#--Covariates z-score scaling
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}


#--Model fit metrics fucntion ----------------------------
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#Mean Absolute Error
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  IMPRECISION = sd(residual, na.rm=T)
  In_IC = mean(obs<upper & obs> lower, na.rm=T)*100
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 IMPRECISION = IMPRECISION,
                 Accuracy = In_IC)
  return(output)
}


#####------------------VARIABLES RECODING/ MORE EXPLORATIONS-----------------------------
###-------------------------------------------------------------------------------
names(covDF2)
datt <- covDF2 #-------rename data
mod_covs <- 34:85 ; length(34:85)#----COVARIATES COLUMNS in the dataset

head(datt[, mod_covs]) #-----view the covariates columns to confirm
datt[,mod_covs] <- apply(datt[,mod_covs], 2, stdize)   #---apply the z-score scaling
head(datt[, mod_covs])   #-----view to confirm
head(datt); names(datt)

#--Explore Observed Count and Building intensity
par(mfrow=c(1,1))
hist(datt$POPN)  #---observed pop total per CU
boxplot(datt$POPN) 
hist(datt$BLDG21) #---Building intensity
plot(datt$POPN, datt$BLDG21) #--scatter
sum(datt$POPN, na.rm=T) #----observed total = 7050773
#------------------------------------------------------------------------
colnames(datt)[mod_covs] <- paste0("x", 1:52, sep="") #--rename covariates 
names(datt)  #------confirm



#----Divide into training and prediction set for covariates selection
dim(dat_train <- datt[!is.na(datt$POPN),]) #-- Training set - all CUs with Observations only (16872 CUs)
dim(dat_pred <- datt[is.na(datt$POPN),])   #-- Prediction set - all unsampled CUs (15095 CUs)



#-----------------------------------------------------------------------------
#---------------PREPARE DATA FOR INLA-SPDE------------------------
#----------------------------------------------------------------
#----Builf the MESH-------------------------
coords = cbind(datt$lon, datt$lat)
coords_train = cbind(dat_train$lon, dat_train$lat)
coords_pred = cbind(dat_pred$lon, dat_pred$lat)

##-View the study locations
plot(coords, cex=0.1, col="black")
plot(coords_train, cex=0.1, col="blue",
     ylab="Latitude",
     xlab="Longitude")
points(coords_pred, cex = 0.1, col="red")
legend("topright", legend = c("Observed", "No Data"),
       col=c("blue", "red"),
       cex=c(0.1,0.1), bty="n")


###----GGPLOT of Observed and missing data locations---------------------------------
plot.dat1 <- data.frame(coords_train)
plot.dat2 <- data.frame(coords_pred)
plot.dat1$status <- rep("Observed", nrow(plot.dat1))
plot.dat2$status <- rep("No Data", nrow(plot.dat2))
plot.dat <- rbind(plot.dat1, plot.dat2); names(plot.dat)
names(plot.dat) <- c("Longitude", "Latitude", "Status")
names(plot.dat)

##----
data_locs <- ggplot(plot.dat, aes(x = Longitude, y = Latitude, color=Status)) +
  annotation_map(map_data("world"), colour = "light green", fill = "dark grey") +
  geom_point(size=2) +
  labs(color = "Status")+
  scale_color_manual(values=c('light green', 'dark red'))+
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
        legend.text = element_text(size = 16, color = "black"),
        legend.title = element_text(size = 16, face = "bold.italic"),
        text=element_text(size=16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)
  )

ggsave("PNG survey locations.png", data_locs, path = paste0(results_path,"/"))



###----Build non-hull mesh
bnd <- inla.nonconvex.hull(as.matrix(coords),-0.03, -0.05, resolution=c(100,100))
mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.6,4), cutoff=0.4)

##---Plot mesh showing study locations---------
par(mfrow=c(1,1))
png(paste0(results_path, "/plots/mesh.png"))
plot(mesh)  
points(coords_train, cex=0.5, col="blue",
     ylab="Latitude",
     xlab="Longitude")
points(coords_pred, cex = 0.5, col="red")
dev.off()
mesh$n #----number of mesh nodes

##---Clear Workspace
rm(dat_pred, coords_pred)

#--------------------------------------------------------------
###---Build projection matrix A
A<-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords));dim(A)

##---Create the SPDE Object
spde <- inla.spde2.matern(mesh, alpha=2)

##----specify the observation indices for estimation 
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)


##---More variables recode
datt$prov <- datt$Prov_ID  #--province
datt$set_typ <- as.factor(as.numeric(datt$TYPE)) #--settlement type

#---settlement - type and province nesting
Zsp <- as(model.matrix( ~ 0 + prov:set_typ, data = datt), "Matrix") 
datt$IDsp <- 1:nrow(datt)
datt$set_prov <- as.factor(apply(Zsp, 1, function(x){names(x)[x == 1]}))#--nesting



##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#---THE 2-STAGE Gamma-Gaussian model---------------------------------------------------------
# In this section, the Gamma-Section model is specified and fitted
# Gaussian Model was fit for the Building Intensity
# While Gamma Model was fit for the Population density
#-------------------------------------------------------------------------------

#------------------------------------------------------------------------------
##--STAGE 1: BUILDING INTENSITY MODEL TO PREDICT FOR CANOPY-COVERED LOCATIONS AS WELL
###-----Building Stack for building intensity---------------------------------
#---specify model covariates and random effects
cov_bld <- datt[,c("x3","x15","x27","x30", "x31","x33","x34", "x35","x36",
                   "x38","x46","x50","x52","set_prov", "set_typ", 
                   "prov", "IDsp")]; dim(cov_bld)

#----stack 
datt$BLDG21 <- datt$BLDG21 + 1# transformed sto have at least a value of 1 
stk_bld <- inla.stack(data=list(y=log(datt$BLDG21)), #the log-transformed response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(cov_bld)
                      ), 
                      #this is a quick name so you can call upon easily
                      tag='est_bld')


#----SPECIFY AND FIT VARIOUS MODELS--------------
#---BEST FIT MODEL for building intensity with decomposed spatial random effects
fbld3<- y ~ -1 + Intercept +  x3 + x15 + x27 + x30 + x31 + x33 + x34 + x35 + x36 + 
  x38 + x46 + x50 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') + f(set_prov, model="iid")

bmod3<-inla(fbld3, #the formula
            data=inla.stack.data(stk_bld,spde=spde),  #the data stack
            family= 'gaussian',   #which family the data comes from
            control.predictor=list(A=inla.stack.A(stk_bld),compute=TRUE),  #compute gives you the marginals of the linear predictor
            control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
            verbose = FALSE) #can include verbose=TRUE to see the log of the model runs

summary(bmod3) #----model summary
bmod3$summary.fix #--extract fixed effects 
bmod3$summary.random #--extract random effects 
bind3 <-inla.stack.index(stk_bld, "est_bld")$data #--estimation indices 
bfit3 <- exp(bmod3$summary.linear.predictor[bind3,"mean"]) #--extract the back transformed predicted response
sum(bfit3)
bfit3L <- exp(bmod3$summary.linear.predictor[bind3,"0.025quant"]) #--extract the back transformed lower bound 
bfit3U <- exp(bmod3$summary.linear.predictor[bind3,"0.975quant"]) #--extract the back transformed upper bound

#--explore the posterior
bld3<- cbind(datt$BLDG21, bfit3)
apply(bld3, 2, sum, na.rm=T)
plot(datt$BLDG21, bfit3, col=c(1,2)) #--scatter plot of observed vs predicted intensity
abline(a=0, b=1)
cor(datt$BLDG21, bfit3)
datt$bld_pred <- bfit3 #--add the predicted intensity to the data


#---Extract and save the beta values----------------------------------
betab <- bmod3$summary.fix #--extract fixed effects 
write.csv(round(betab,5), paste0(results_path,"/updated2/Betas_for_building_intensity_model.csv"))
##--Extract the DIC and WAIC 
mod.fitb<- bmod3$dic$dic


#----------------------------------------------------------------------------
##---STAGE 2: DENSITY MODEL WITH PREDICTED BUILDING INYENSITY
#------------------------------------------------------------------------------

datt$dens2 <- datt$POPN/datt$bld_pred #---DEFINE THE DENSITY VARIABLE
hist(log(datt$dens2), col="brown") #--view the histogram 
datt$dens2[datt$dens2==0] = 0.000001 #--There is at least one individual in the CU's with zero data


####---Density Stack
covars_dens <- datt[,c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
                       "x39","x45","x46","x48","x51","x52","set_prov", "set_typ", 
                       "prov", "IDsp")]; dim(covars_dens)

#---Build the stack for the training set
stk_dens <- inla.stack(data=list(y=datt$dens2), #the response
                       
                       A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                       
                       effects=list(c(list(Intercept=1), #the Intercept
                                      iset),  #the spatial index
                                    #the covariates
                                    list(covars_dens)
                       ), 
                       #this is a quick name so you can call upon easily
                       tag='est_dens')


#----Specify and fit the best model for density------------------------------------
##==========================
form1c <- y ~ -1 + Intercept +  x3 + x11 + x15 + x29 + x30 + x31 + x32 + x35 + x36 + 
  x39 + x45+ x46 + x48 + x51 + x52 + f(spatial.field, model=spde) + f(IDsp, model='iid') +
  f(set_typ, model='iid')

mod1c <-inla(form1c, #the formula
             data=inla.stack.data(stk_dens,spde=spde),  #the data stack
             family= 'gamma',   #which family the data comes from
             control.predictor=list(A=inla.stack.A(stk_dens),compute=TRUE),  #compute gives you the marginals of the linear predictor
             control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
             verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
#summary(mod1c) #----model summary


#---Extract and save the beta values 
betagg <- mod1c$summary.fix
write.csv(round(betagg,5), paste0(results_path,"/updated2/Betas_for_density_model.csv"))


#--extract random effects 
ind1 <-inla.stack.index(stk_dens, "est_dens")$data #--estimation indices 
fit1c <- exp(mod1c$summary.linear.predictor[ind1,"mean"]) #--extract the back transformed density
fit1cU <- exp(mod1c$summary.linear.predictor[ind1,"0.975quant"]) #--extract the back transformed lower bound
fit1cL <- exp(mod1c$summary.linear.predictor[ind1,"0.025quant"]) #--extract the back transformed upper bound

fit11c <- fit1c*datt$bld_pred #----Back transformed pop_hat
fit11cU <- fit1cU*datt$bld_pred #---lower bound
fit11cL <- fit1cL*datt$bld_pred  #----Upper bound 
(sum_fit <- sum(fit1c*datt$bld_pred))
(POPb <- cbind(datt$POPN, fit11c))
apply(POPb, 2, sum, na.rm=T)
plot(datt$POPN, fit11c, col=c("red", "blue")) #--scatter plot
abline(a=0, b=1)
cor(datt$POPN, fit11c) 


##---scatter plot
plot(fit11c,datt$POPN,  xlab = "Observed", 
     ylab = "Predicted", col=c('dark green','orange'),
     pch=c(16,16), cex.axis=1.5)
abline(0,1)
legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
       bty="n", cex=1.5)


##---Extract DIC
(mod.fit<- mod1c$dic$dic)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####-----POSTERIOR SIMULATION-----------------------------------------------------
#-------------------------------------------------------------------------------------

#---Extract Settlement type random effects
#----Function for extracting settlement type random effects
set_t <- function(dat, st)
{
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$set_typ2[i] = st[j]
    }
    
  }
  dat$set_typ2
}

datt$set_typ2 <- set_t(datt, st=mod1c$summary.random$set_typ$mean)
#class(datt$set_typ2)

#----Posterior Simularion-------------------------------------------
simPops_gg <- function(model, dat, Aprediction, run)
{
  fixedeff  <- dens_hat <- pop_hat <- matrix(0, nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =  481561959 #--simulation seed to ensure reproducible estimates
  set.seed(inla.seed)
  print(inla.seed)
  m1.samp <- inla.posterior.sample(run, model, seed = inla.seed ,selection=list(x3=1, x11=1, x15=1,
                                                                                x29=1, x30=1, x31=1,
                                                                                x32=1, x35=1, x36=1,
                                                                                x39=1, x45=1, x46=1,
                                                                                x48=1, x51=1, x52=1),num.threads="1:1")
  
  sfield_nodes_mean <- model$summary.random$spatial.field['mean']
  field_mean <- (Aprediction%*% as.data.frame(sfield_nodes_mean)[, 1])
  for(i in 1:run)
  {
    fixedeff[,i] <-  
      model$summary.fixed['Intercept', 'mean'] +
      m1.samp[[i]]$latent[1,] * dat[,'x3'] +
      m1.samp[[i]]$latent[2,] * dat[,'x11'] +
      m1.samp[[i]]$latent[3,] * dat[,'x15'] +
      m1.samp[[i]]$latent[4,] * dat[,'x29'] +
      m1.samp[[i]]$latent[5,] * dat[,'x30'] + 
      m1.samp[[i]]$latent[6,] * dat[,'x31'] +
      m1.samp[[i]]$latent[7,] * dat[,'x32'] +
      m1.samp[[i]]$latent[8,] * dat[,'x35'] +
      m1.samp[[i]]$latent[9,] * dat[,'x36'] +
      m1.samp[[i]]$latent[10,] * dat[,'x39'] +
      m1.samp[[i]]$latent[11,] * dat[,'x45'] +
      m1.samp[[i]]$latent[12,] * dat[,'x46'] +
      m1.samp[[i]]$latent[13,] * dat[,'x48'] +
      m1.samp[[i]]$latent[14,] * dat[,'x51'] +
      m1.samp[[i]]$latent[15,] * dat[,'x52'] +
      model$summary.random$IDsp['mean'][,1] +
      dat$set_typ2 +
      rnorm(nrow(dat), 0, 1/m1.samp[[i]]$hyperpar[1]) +
      
      field_mean[,1]
    
    dens_hat[,i]<- exp(fixedeff[,i])
    pop_hat[,i] <- dens_hat[,i]*dat$bld_pred
  }
  
  mean_dens_hat <- apply(dens_hat, 1, mean, na.rm=T) #
  mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
  median_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.5), na.rm=T) #
  lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  uncert_pop_hat <- (upper_pop_hat - lower_pop_hat)/mean_pop_hat#
  
  
  dat$mean_dens_hat <- mean_dens_hat
  dat$mean_pop_hat <- mean_pop_hat
  dat$median_pop_hat <- median_pop_hat
  dat$lower_pop_hat <- lower_pop_hat
  dat$upper_pop_hat <- upper_pop_hat
  dat$uncert_pop_hat <- uncert_pop_hat
  dat$sd_pop_hat <- sd_pop_hat
  
  
  
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
  
}

#rm(dat.sim, sim.pops_nb, sim.pops, llg.est, prov.est2)
run=2000 #--number of simulations
system.time(str(sim.pops_gg <- simPops_gg(mod1c,datt, A, run)))



#---------------------------------------------------------
####---View a random subset of the Posterior samples
#---------------------------------------------------------

##-----Traceplots for posterior samples
samp <- sample(nrow(datt), 6)
par(mfrow=c(3,2), mar=c(5,5,2,1))
for(j in samp)
{
  plot.ts(sim.pops_gg$pop_hat[j,], ylab = "pop_hat",
          cex.main=2,  main=paste0("Pop_hat samples for CU ", datt$CU_Name[j], sep=""))
  abline(h=mean(sim.pops_gg$pop_hat[j,]), lwd=2, col=2)
}


##---Histograms of posterior samples
for(j in samp)
{
  hist(sim.pops_gg$pop_hat[j,],xlab="pop_hat", prob=T, ylim=c(0,0.01),
       cex.axis=1.5, main=paste0("Pop_hat histogram for CU ", 
                                 datt$CU_Name[j], sep=""))
  abline(v=mean(sim.pops_gg$pop_hat[j,]), lwd=2, col="red")
}




#####----Join the posterior draws to the dataframe------------------------------------------------------
sum(datt$pop_hat1 <- fit11c)
dim(datt)

data.all <- datt[,c("Prov_Name" ,"Dist_Name", "LLG_Name", "Ward_Name", "pop_hat1")]
data.all <- data.frame(cbind(data.all, sim.pops_gg$pop_hat))#---Contains all the simulated posterior matrix
dim(data.all);names(data.all)


####-----------ADMIN TOTALS AND UNCERTAINTIES----------------------------------------
thin=5 #----Thinning chooses every 5th sample for posterior inference 
thinn <- seq((6+0.2*run),run, by=thin) #--There is a burn-in period of the first 20% of the total sample


#----NATIONAL TOTAL AND UNCERTAINTIES------------------------------
nat_total <- function(dat, thinn)
{
  p_hat <- dat[,thinn]
  tots <- apply(p_hat,2, sum, na.rm=T) #Col sums
  
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean, lower=tot_lower, median=tot_median, upper=tot_upper))))
}
(national_gg <- nat_total(data.all, thinn))

#--Save data
write.csv(national_gg, file=paste0(results_path, "/updated2/National_estimates_main_gamma-gaussian.csv"))


##----PROVINCIAL TOTALS AND UNCERTAINTIES----------------------------------------
prov_est_gg <- function(datp, thinn)
{
  provnames <- unique(datp$Prov_Name)
  outP <- matrix(0, nrow=length(provnames), ncol=4)
  for(j in 1:length(provnames))
  {
    prov <- datp[datp$Prov_Name==provnames[j],]
    ptots <- apply(prov[,thinn], 2, sum, na.rm=T)
    ptots_sd <- sd(ptots, na.rm=T)
    
    ptot_mean1  <- mean(ptots, na.rm=T)
    ptot_lower <- quantile(ptots, probs=c(0.025))
    ptot_meadian <- quantile(ptots, probs=c(0.5))
    ptot_upper <- quantile(ptots, probs=c(0.975))
    
    pestimates <- round(c(ptot_mean1, ptot_lower, ptot_meadian, ptot_upper), 2)
    outP[j,] <- pestimates
  }
  outP <- data.frame(outP)
  return(PROV_est <- data.frame(names = provnames,
                                total = outP[,1],
                                lower = outP[,2],
                                median = outP[,3],
                                upper = outP[,4]))
}
(prov.est_gg <- prov_est_gg(data.all,thinn))
sum(prov.est_gg$total)

#--save data 
write.csv(prov.est_gg, file=paste0(results_path, "/updated2/Province_estimates_main_gamma-gaussian.csv"))



###---DISTRICT LEVEL TOTALS AND UNCERTAINTIES----------------------------------------
dist_est_gg <- function(datd, thinn)
{
  dnames <- unique(datd$Dist_Name)
  outd <- matrix(0, nrow=length(dnames), ncol=4)
  
  for(j in 1:length(dnames))
  {
    dist <- datd[datd$Dist_Name==dnames[j],]
    dtots <- apply(dist[,thinn], 2, sum, na.rm=T)
    
    dtots_sd <- sd(dtots, na.rm=T)
    
    dtot_mean1  <- mean(dtots, na.rm=T)
    
    dtot_lower <- quantile(dtots, probs=c(0.025))
    dtot_median <- quantile(dtots, probs=c(0.5))
    dtot_upper <- quantile(dtots, probs=c(0.975))
    destimates <- round(c(dtot_mean1, dtot_lower, dtot_median, dtot_upper), 2)
    outd[j,] <- destimates
  }
  outd <- data.frame(outd)
  return(dist_est <- data.frame(names = dnames,
                                total = outd[,1],
                                lower = outd[,2],
                                median = outd[,3],
                                upper = outd[,4]))
}

(dist.est_gg <- dist_est_gg(data.all,thinn))
sum(dist.est_gg$total)

#--Save data
write.csv(dist.est_gg, file=paste0(results_path, "/updated2/Districts_estimates_main_gamma-gaussian.csv"))


###---LLG level TOTALS AND UNCERTAINTIES ------------------------### 
#----------------------------------------------------------
llg_est_gg <- function(datl, thinn)
{
  lnames <- unique(datl$LLG_Name)
  outl <- matrix(0, nrow=length(lnames), ncol=4)
  
  for(j in 1:length(lnames))
  {
    llg <- datl[datl$LLG_Name==lnames[j],]
    ltots <- apply(llg[,thinn], 2, sum, na.rm=T)
    ltots_sd <- sd(ltots, na.rm=T)
    
    ltot_mean1  <- mean(ltots, na.rm=T)
    ltot_lower <- quantile(ltots, probs=c(0.025))
    ltot_median <- quantile(ltots, probs=c(0.5))
    ltot_upper <- quantile(ltots, probs=c(0.975))
    
    lestimates <- round(c(ltot_mean1, ltot_lower, ltot_median, ltot_upper), 3)
    outl[j,] <- lestimates
  }
  outl <- data.frame(outl)
  return(llg_est <- data.frame(names = lnames,
                               total = outl[,1],
                               lower = outl[,2],
                               median = outl[,3],
                               upper = outl[,4]))
}

(llg.est_gg <- llg_est_gg(data.all,  thinn))
sum(llg.est_gg$total)

#--Save data
write.csv(llg.est_gg, file=paste0(results_path, "/updated2/LLG_estimates_main_gamma_gaussian.csv"))


###---CU level TOTALS AND UNCERTAINTIES------------------------------
#----------------------------------------------------------
cu_est_gg <- function(dat, dat2, thinn)
{
  dat$mean  <- apply(dat2[,thinn], 1, mean, na.rm=T)
  dat$lower <- apply(dat2[,thinn], 1, quantile, probs=c(0.025), na.rm=T)
  dat$median <- apply(dat2[,thinn], 1, quantile, probs=c(0.5), na.rm=T)
  dat$upper <- apply(dat2[,thinn], 1, quantile, probs=c(0.975), na.rm=T)
  dat$uncertainty <- (dat$upper-dat$lower)/dat$mean
  return(dat)
}
cu.est <- cu_est_gg(sim.pops_gg$est_data, data.all,thinn)
head(cu.est)
sum(cu.est$mean, na.rm=T)
Vars2Include <- c("CLUSTER_ID", "Uniq_ID", "Prov_Name" ,"Dist_Name", "LLG_Name", "Ward_Name",
                  "CU_Name", "TYPE", "lon", "lat", "SOURCE", "POPN", "BLDG21", "mean", "median",
                  "lower", "upper", "uncertainty")
dim(cu.data <- cu.est[,Vars2Include]); head(cu.data)

#---save data
write.csv(cu.data, file=paste0(results_path, "/updated2/CU_data_with_estimates.csv"))
#---------------------------------------------------------------------------------------------
#--END OF  MODEL FITTING, POSTERIOR SIMULATION AND POSTERIOR EVALUATION-----------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#++++OPTIONAL++++++++++++++++++++++++++++++++++++++++++++++
####-----CROSS-VALIDATION -----------------------------------------

#--Load the Cross Validation script 
source(paste0(data_path, "gg_model_scripts/gg_model_cross_validation.R"))
#print(nm)
k_folds=5  #-5 folds cross validation
(cv <- cross_val(datt, mod1c, fit11c, mesh, spde,
                 A, shp, form1c, k_folds))

cv$metrics #---extract model fit metrics

str(cv$data)


cv.data <- cv$data
dim(cv_data <-  do.call(rbind, cv.data)) # join the matrcices in the list by row as a single matrix 
names(cv_data)

table(cv_data$fold)


dim(cv11 <- cv_data[cv_data$fold==1, ])

plot(cv11$pred, cv11$test)
abline(a=0, b=1)
#save the model fit metrics
write.csv(cv$metrics, paste0(results_path, "/cross_validation.csv"))




###-----------MORE POSTERIOR EXPLORATIONS AND PLOTS -------------------------------
#----Add to CU shapefile
shp2 <- shp
shp2$mean <- cu.data$mean
shp2$lower <- cu.data$lower
shp2$median <- cu.data$median
shp2$upper <- cu.data$upper
shp2$uncertainty <- cu.data$uncertainty


#require(mapview) #interactive maps

mapview(shp2, zcol="mean")


##----Save workspace
save.image(paste0(results_path, "/workspace_gg_model.Rdata"))