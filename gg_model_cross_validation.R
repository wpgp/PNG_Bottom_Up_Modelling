#

####--TITLE: R SCRIPTS FOR MODELLED POPULATION ESTIMATES FOR PAPUA NEW GUNEA BASED ON R-INLA----#####
####--METHODS: GEOSTATISTICAL BAYESIAN HIERARCHICAL REGRESSION MODEL
####--AUTHOR: DR CHIBUZOR CHRISTOPHER NNANATU
####--INSTITUTION: WORLDPOP, UNIVERSITY OF SOUTHAMPTON 
####--DATE:19 DECEMBER 2022.
#nv <- sample(1000000,1)
#print(nv)
par(mfrow=c(5,2), mar=c(2,2,1,1))
cross_val <- function(dat, mod, mod_pred, mesh, spde,
                      A, shp, formula, k_folds)
{
  set.seed(536444)   
  N <- nrow(dat)
  ######
  ind_train <- factor(sample(x = rep(1:k_folds, each = floor((N+1)/ k_folds)),  # Sample IDs for training data
                             size = N))
  
  table(as.numeric(ind_train)) 
  dat$k_fold <- as.numeric(ind_train)
  coords <- cbind(dat$lon, dat$lat)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  measures <- c("MAE", "RMSE", "Abias","IMPRECISION","%Accuracy")
  metrics_cv <- matrix(0, ncol=length(k_uniq), nrow=length(measures))#
  

 # plot(mod_pred, dat$POPN, xlab = "Predicted", 
      # ylab = "Observed", col=c('dark green','orange'),
      # pch=c(16,16), cex.axis=1.5)
 # abline(0,1)
  #legend("topleft", c("Predicted", "Observed"), col=c("dark green", "orange"), pch=c(16,16),
        # bty="n", cex=1.5)
  
  cvs <- list()
  for(i in 1:length(k_uniq))
  {
    k_uniq <-sort(unique(dat$k_fold))
    
    print(paste0("fold_", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[-train_ind, ]) #---test set for fold i
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[-train_ind,]
    
    
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    png(paste0(results_path,"/fold_",i,"_sample points for cross-validation.png"))
    plot(shp)
    points(train_coords, pch=16,  col="blue", cex=0.6)
    points(test_coords, pch=2,  col="red", cex=0.6)
    dev.off()
    
    
    ########################
    covars_train <- train[,c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
                             "x39","x45","x46","x48","x51","x52","set_prov", "set_typ", 
                             "prov", "IDsp")]; dim(covars_train)
    
    #---Build the stack for the training set
    stk_train <- inla.stack(data=list(y=train$dens2), #the response
                            
                            A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars_train)
                            ), 
                            #this is a quick name so you can call upon easily
                            tag='train')
    
    
    
    covars_test <- test[,c("x3","x11","x15","x29","x30", "x31","x32","x35", "x36",
                           "x39","x45","x46","x48","x51","x52","set_prov", "set_typ", 
                           "prov", "IDsp")]; dim(covars_test)
    ###################################
    
    ##==========================
    #formula = form1c
    ###---Rerun INLA for model test prediction
    model <-inla(formula, #the formula
                 data=inla.stack.data(stk_train,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model)
    
    
    
    #mod = mod1c
    #######################----Extract Spatial Random effects
    sfield_nodes_mean <- mod$summary.random$spatial.field['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
    sfield_nodesL <- mod$summary.random$spatial.field['0.025quant']
    fieldL <- (A%*% as.data.frame(sfield_nodesL)[, 1])
    
    sfield_nodesU<- mod$summary.random$spatial.field['0.975quant']
    fieldU <- (A%*% as.data.frame(sfield_nodesU)[, 1])
    
    
    
    ###----Extract settlement type random effects
    set <- set_t(dat, mod$summary.random$set_typ$mean)
    setL <- set_t(dat, mod$summary.random$set_typ$`0.025quant`)
    setU <- set_t(dat, mod$summary.random$set_typ$`0.975quant`)
    
    
    ##--------
    fixed <-  
      model$summary.fixed['Intercept', 'mean'] +
      model$summary.fixed['x3', 'mean']* test[,'x3'] +
      model$summary.fixed['x11', 'mean'] * test[,'x11'] +
      model$summary.fixed['x15', 'mean'] * test[,'x15'] +
      model$summary.fixed['x29', 'mean'] * test[,'x29'] +
      model$summary.fixed['x30', 'mean'] * test[,'x30'] + 
      model$summary.fixed['x31', 'mean'] * test[,'x31'] +
      model$summary.fixed['x32', 'mean'] * test[,'x32'] +
      model$summary.fixed['x35', 'mean'] * test[,'x35'] +
      model$summary.fixed['x36', 'mean'] * test[,'x36'] +
      model$summary.fixed['x39', 'mean'] * test[,'x39'] +
      model$summary.fixed['x45', 'mean'] * test[,'x45'] +
      model$summary.fixed['x46', 'mean'] * test[,'x46'] +
      model$summary.fixed['x48', 'mean',] * test[,'x48'] +
      model$summary.fixed['x51', 'mean'] * test[,'x51'] +
      model$summary.fixed['x52', 'mean'] * test[,'x52'] +
      
      mod$summary.random$IDsp['mean'][-train_ind,1] +
      set[-train_ind] +
      rnorm(nrow(test),0, 1/model$summary.hyperpar$mean[1]) + #
      field_mean[-train_ind,1]
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld_pred)
    
    
    
    
    ######----Lower
    fixedL <-  
      model$summary.fixed['Intercept', '0.025quant'] +
      model$summary.fixed['x3', '0.025quant']* test[,'x3'] +
      model$summary.fixed['x11', '0.025quant'] * test[,'x11'] +
      model$summary.fixed['x15', '0.025quant'] * test[,'x15'] +
      model$summary.fixed['x29', '0.025quant'] * test[,'x29'] +
      model$summary.fixed['x30', '0.025quant'] * test[,'x30'] + 
      model$summary.fixed['x31', '0.025quant'] * test[,'x31'] +
      model$summary.fixed['x32', '0.025quant'] * test[,'x32'] +
      model$summary.fixed['x35', '0.025quant'] * test[,'x35'] +
      model$summary.fixed['x36', '0.025quant'] * test[,'x36'] +
      model$summary.fixed['x39', '0.025quant'] * test[,'x39'] +
      model$summary.fixed['x45', '0.025quant'] * test[,'x45'] +
      model$summary.fixed['x46', '0.025quant'] * test[,'x46'] +
      model$summary.fixed['x48', '0.025quant',] * test[,'x48'] +
      model$summary.fixed['x51', '0.025quant'] * test[,'x51'] +
      model$summary.fixed['x52', '0.025quant'] * test[,'x52'] +
      
      mod$summary.random$IDsp['0.025quant'][-train_ind,1] +
      setL[-train_ind] + 
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.025quant`[1]) + #
      fieldL[-train_ind,1]
    
    dens_htL <- exp(fixedL)
    sum(pop_htL <- dens_htL*test$bld_pred)
    
    
    #=========Upper
    fixedU <-  
      model$summary.fixed['Intercept', '0.975quant'] +
      model$summary.fixed['x3', '0.975quant']* test[,'x3'] +
      model$summary.fixed['x11', '0.975quant'] * test[,'x11'] +
      model$summary.fixed['x15', '0.975quant'] * test[,'x15'] +
      model$summary.fixed['x29', '0.975quant'] * test[,'x29'] +
      model$summary.fixed['x30', '0.975quant'] * test[,'x30'] + 
      model$summary.fixed['x31', '0.975quant'] * test[,'x31'] +
      model$summary.fixed['x32', '0.975quant'] * test[,'x32'] +
      model$summary.fixed['x35', '0.975quant'] * test[,'x35'] +
      model$summary.fixed['x36', '0.975quant'] * test[,'x36'] +
      model$summary.fixed['x39', '0.975quant'] * test[,'x39'] +
      model$summary.fixed['x45', '0.975quant'] * test[,'x45'] +
      model$summary.fixed['x46', '0.975quant'] * test[,'x46'] +
      model$summary.fixed['x48', '0.975quant',] * test[,'x48'] +
      model$summary.fixed['x51', '0.975quant'] * test[,'x51'] +
      model$summary.fixed['x52', '0.975quant'] * test[,'x52'] +
      
      mod$summary.random$IDsp['0.975quant'][-train_ind,1] +
      setU[-train_ind] + 
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$`0.975quant`[1]) + #
      fieldU[-train_ind,1]
    
    dens_htU <- exp(fixedU)
    sum(pop_htU <- dens_htU*test$bld_pred)
    
    
    ###
    
    plot(test$POPN, pop_ht, xlab = "Observed", 
         ylab = "Predicted", col=c('dark green','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("dark green", "orange"), pch=c(16,16),
           bty="n", cex=1.5) 
    
    
    (met <- model_metrics(test$POPN,  
                          pop_ht, pop_htU,  pop_htL))
    metrics_cv[,i] <- as.vector(unlist(met)) 
    
    cvs[[i]] <- data.frame(test=test$POPN, pred=pop_ht, fold=rep(k_uniq[i], length(as.vector(pop_ht))))
  }
  
  
  stat = data.frame(metrics =c("MAE", "RMSE", "Abias","IMPRECISION","%Accuracy"))
  values = data.frame(metrics_cv)
  colnames(values) = paste0("Fold_", 1:k_folds, sep="")
  metrics = bind_cols(stat, values)
  
  return(list(metrics= metrics, data = cvs))
}

#print(nm)
#k_folds=10

#(cv <- cross_val(datt, mod1c, fit11c, mesh, spde,
 #                A, shp, form1c, k_folds))





#cv$metrics
#write.csv(cv, paste0(results_path, "/cross_validation.csv"))



#mean(c(11781559, 11998136, 11852921, 11745906))

# Basic scatter plot.
#dt1 <- data.frame(test = datt$POPN, pred=fit11c, fold=rep(0, nrow(datt)))
#dtt <- bind_rows(dt1, cv$data[[1]], cv$data[[2]],cv$data[[3]],
 #                cv$data[[4]],cv$data[[5]])

#dtt$fold1 <- factor(dtt$fold)
#levels(dtt$fold1) <- paste0("fold", 0:5, sep="")
#table(dtt$fold1)

# Library

#library(hrbrthemes)

#p1 <- ggplot(dtt, aes(x=test, y=pred)) + 
  #geom_point( color=c("69b3a2"), size=3) +
  #geom_smooth(method=lm , color="orange", fill=c("#69b3a2"), se=TRUE) +
  #geom_smooth(method=lm , color="red", se=FALSE) +
  #labs(
   # x = "Observed Count (Test)",
  #  y = "Predicted Count"
  #) +
  #theme(#panel.background = element_rect(fill = "#006994", colour = "black"),
  #  legend.position = "top", legend.spacing.x = unit(1.0, "cm"),
  #  legend.text = element_text(size = 18, color = "black"),
  #  legend.title = element_text(size = 20, face = "bold.italic"),
  #  text=element_text(size=24),
   # axis.title.x = element_text(size = 20),
  #  axis.title.y = element_text(size = 20)
 # )+
 # ggtitle("Out of Sample Cross-Validation")+
 # facet_wrap(dtt$fold1)

#p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

