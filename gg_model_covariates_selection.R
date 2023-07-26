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

#-------------------------------------------------------------------------------
#------------------------Load datasets-----------------------------------------        
#--------------------------------------------------------------------------------
shp <- readOGR(dsn = paste0(data_path,"gg_model_input_data/Boundary_files/CU"), layer = "PNG_CU_32100_B") #combined shapefile 
#plot(shp);dim(shp); names(shp)

##load Rdata--
load(paste0(data_path,"gg_model_input_data/Survey_data/ModelDataPNG_BI.RData")) #-
#-----------------------------------------------------------------------------
#--------------------------------------------------------------------------
###---Explore  and clean data ----------------------------------------
#--------------------------------------------------------------------------
ls()
names(covs)
covDF <- as.data.frame(covs)
#plot(covDF$lon, covDF$lat)
#names(covDF); str(covDF); dim(covDF); dim(shp)# shp has one extra row
shp <- shp[-32100,] #--remove row 32100 from the shapefile as it does not exist
dim(covDF); dim(shp) #---confirm
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
##----Extract the coordiates - the centroids of the CUs
#--------------------------------------------------------------------
library(sp)
covDF$lon <- coordinates(shp)[,1]#--add lon-lat to the data
covDF$lat <- coordinates(shp)[,2]
covDF2 <- covDF
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#----------------------Covariates scaling (Z-score) function--------
#----------------------------------------------------------------------------
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}
#-----------------------------------------------------------------------------


#---------------- Model fit metrics fucntion ----------------------------
model_metrics <- function(obs, pred, upper, lower)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  IMPRECISION = sd(residual, na.rm=T)
  In_IC = mean(obs<upper & obs> lower)*100
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 BIAS = abs(BIAS),
                 IMPRECISION = IMPRECISION,
                 Accuracy = In_IC)
  return(output)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####------------------VARIABLES RECODING/EXPLORATION-----------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
names(covDF2)
datt <- covDF2 #-------rename data
mod_covs <- 34:85 ; length(34:85)#----COVARIATES COLUMNS in datt

head(datt[, mod_covs]) #-----view the covariates columns
datt[,mod_covs] <- apply(datt[,mod_covs], 2, stdize)   #---apply the z-score scaling
head(datt[, mod_covs])   #-----view to confirm
head(datt); names(datt)

#-------
par(mfrow=c(1,1))
hist(datt$POPN)  #---observed pop total per CU
boxplot(datt$POPN) 
hist(datt$BLDG21) #---Building intensity
plot(datt$POPN, datt$BLDG21) #--scatter
sum(datt$POPN, na.rm=T) #----observed total = 7050773
#------------------------------------------------------------------------
colnames(datt)[mod_covs] <- paste0("x", 1:52, sep="") #--rename covariates 
names(datt)  #------confirm

dim(dat_train <- datt[!is.na(datt$POPN),]) #-- Training set - all CUs with Observations only (16872 CUs)
dim(dat_pred <- datt[is.na(datt$POPN),])   #-- Prediction set - all unsampled CUs (15095 CUs)

#----------------Covariates Selection-----------------------------------
#install.packages("car")
library(car) ##--For calculating variance inflation factor (vif)
library(dplyr)
library(tidyverse)

#######-------For BLD INTENSITY================================================================
dat_train$lbld <- log(dat_train$BLDG21)
dim(covs_bld <- dat_train[,c(91, 95,mod_covs)])#--subset for variables selection
covs_bld1 <- covs_bld %>% drop_na() #- model covariates only without NAs
dim(covs_bld2 <- covs_bld1 %>% dplyr::select(!POPN))
dim(covs_bld2 <- covs_bld2[is.finite(covs_bld2$lbld),]) #--checks

#-----------------FITTIG THE STEPWISE REG MODELS-------------------------------------
fbld <- glm(lbld ~., data=covs_bld2, family = gaussian(link=identity)) #---Gaussian GLM

step_bld <- stepAIC(fbld, scale = 0,
                    direction = c("both"),#--SELECTION OF COVARIATES IN BOTH DIRECTIONS
                    trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                    k = 2)
step_bld


##----BFurther checks on the best covariates for significance and multicollinerity 
##---to avoid overfititng and variance inflation

fbld1 <- glm(formula = lbld ~ x1 + x2 + x3 + x4 + x5 + x6 + x8 + x9 + 
      x10 + x12 + x13 + x14 + x15 + x17 + x18 + x19 + x21 + x22 + 
      x23 + x24 + x26 + x27 + x30 + x31 + x33 + x34 + x35 + x36 + 
      x37 + x38 + x41 + x43 + x46 + x47 + x49 + x50 + x52, family = gaussian(link = identity), 
   data = covs_bld2)

##---Calculate the VIF (multicollinearity)and retain only covariates with VIF < 5
#-------------------------------------------------------------------------------
vif_bld = vif(fbld1)
vif_bld[which(vif_bld < 5)]

##-----Refit model with covariates with VIF < 5 and retain only statistically significant covs
fbld2 <- glm(formula = lbld ~ x3 + x15 + x27 + x30 + x31 + x33 + x34 + x35 + x36 + 
               x38 + x46 + x50 + x52, family = gaussian(link = identity), 
             data = covs_bld2)
summary(fbld2)

#-------------------------------------
#--Final covariates selected
cov_bld <- c("x3","x15","x27","x30", "x31","x33","x34", "x35","x36",
             "x38","x46","x50","x52")

cov_names <- names(covDF)[mod_covs]
covariates_bld <- cov_names[c(3, 15, 27, 30, 31, 33, 34, 
                              35, 36, 38, 46, 50, 52)]
bld_covs <- data.frame(cov = cov_bld, name=covariates_bld)
#-------------------------------------

##----Make and save correlation plot
require(corrplot)
png(paste0(results_path, "/cor_plots for building intensity.png"))
corrplot(
  cor(covs_bld[,cov_bld]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
dev.off()


##--------------------------------------------------------------------
##-----COVARIATES SELECTION FOR DENSITY-----------------------------------------------------
##-----------------------------------------------------------------------
dim(covs_dens <- covs_bld1[is.finite(covs_bld1$lbld),]) #--select only finite rows 
covs_dens$bld_prd <- round(exp(predict(fbld2))) #- model covariates only without NAs
covs_dens$dens <- covs_dens$POPN/covs_dens$bld_prd #--calculate density
covs_dens$ldens <- log(covs_dens$dens) #--calculate density
dim(covs_dens1 <- covs_dens%>% dplyr::select(!c(lbld, bld_prd, dens)))
dim(covs_dens1 <- covs_dens1[is.finite(covs_dens1$ldens),]) #--select only finite rows 

hist(covs_dens1$ldens, col="brown") #--histogram of density
covs_dens1$ldens[covs_dens1$ldens==0] = 0.000001 #set CU with zero building intensity to have at least  individual per 10^6 intensity
                                            ##4918 16520 16521 16527 18142 18149 18262 18267 18297 18414 20900 20908; 12 CU's no observations


#-----------------FITTIG THE STEPWISE REG MODELS-------------------------------------
fdens <- glm(ldens ~., data=covs_dens1) #---negative binomial

step_dens <- stepAIC(fdens, scale = 0,
                    direction = c("both"),
                    trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                    k = 2)
#step_dens


##--------------------------Refit with covariates with VIF < 5 and retain 


##---NOTE: All coVariates selected for other two competing models - nbinomial and Quasi-Poisson
##---were retained to ensure they could be comapared
#-------------------------------------
cov_dens <- c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
             "x39","x45","x46","x48","x51","x52")                   
                                                                   
cov_names <- names(covDF2)[mod_covs]
covariates_dens <- cov_names[c(3, 11, 15, 29, 30, 31, 32, 
                              35, 36, 39, 45, 46, 48, 51, 52)]
dens_covs <- data.frame(cov = cov_dens, name=covariates_dens)
#-------------------------------------

##----Make correlation plot
require(corrplot)
png(paste0(results_path, "/cor_plots_for_density.png"))
corrplot(
  cor(covs_dens[,cov_dens]),
  method = 'square',
  type = 'upper',
  tl.col = 'black',
  tl.cex = 1,
  col = colorRampPalette(c('purple', 'dark green'))(200)
)
dev.off()
