#####################################################
# Miranda Tao
# created: 8/24/18
# updated: 9/18/18
# Model selection with cross validation on mixed effects models
####################################################
rm(list = ls())

my_lib <- '~/R/my_library/'
#install.packages("glmnet", repos = "http://cran.us.r-project.org",lib = my_lib)
# install.packages("doBy", lib = my_lib)
# install.packages('corrplot',lib = my_lib)
# install.packages('MuMIn', lib = my_lib)
# install.packages('AICcmodavg',lib = my_lib)

pacman::p_load(data.table, dplyr,tidyr,stringr)


library(glmnet,lib.loc = my_lib)
library(doBy,  lib.loc = my_lib)
library(gtools)
library(corrplot,lib.loc =my_lib) # correlation plot
library(lme4)
# library(MuMIn, lib.loc = my_lib) # Multi-Model Inference - MODEL SELECTION
library(AICcmodavg, lib.loc = my_lib) # compute AIC for mixed effect models
source("http://www.sthda.com/upload/rquery_cormat.r")
require(ggplot2)

#===================================================
# data description:
#     - prepped data of malaria health expenditure along with 19 covariates we gathered
#     - The data consists of a subset of sub-saharan african countries, from 2000 to 2015
#     - Data are extracted from multiple online sources and compiled toegther
#     - malaria covariates are pulled from internal database from the Global Burden of Disease study

cov_dt <- fread('~/covariates.csv')
data_ssa <- fread('~/data_ssa.csv')

# merge data with covariates
data_ssa  <- merge(cov_dt, data_ssa, by = c('location_id','year_id'))

# get names for all covariates as predictors
predictors <- names(cov_dt)[3:22]

# plot a correlation matrix
col<- colorRampPalette(c("blue", "white", "red"))(20)
cormat<-rquery.cormat(data_ssa[,predictors,with=F], type="upper", col=col)

#*****************************************************************************************************************
#                                              find the best lambda
#*****************************************************************************************************************

x <- model.matrix(as.formula(paste0('~',paste(predictors,sep = ' ',collapse = '+'))), data_ssa)
y <- data_ssa$y

# find the smallest lambda
cvfit = cv.glmnet(x, y)
lambda_min = cvfit$lambda.min

# run glmnet and zero out coefficients
fit <- glmnet(x,y, family = "gaussian",lambda = lambda_min)
coef(fit)

# variables with zero coefficients need to be dropped
drop <- c('logit_pfpr','logit_itn_pfpr_interaction','logit_anc1','logit_sba','logit_sdi')

#*****************************************************************************************************************
#                                         MODEL SELECTION - with Cross validation
#*****************************************************************************************************************

## Create vectors for outcome and predictors
y    <- c("y")
predictors <- predictors[!(predictors %in% drop)]
mix <- '+(1|location_name)+(1|region_name)'
dataset <- data_ssa

#========================================
## STEP1: Create list of models
#=======================================
list.of.models <- lapply(seq_along((predictors)), function(n) {
  
  LHS <- y
  RHS <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")
  
  paste(LHS, paste0(RHS, mix), sep = "  ~  ")
})

## Convert to a vector
vector.of.models <- unlist(list.of.models)

#==========================================================================
# STEP 2: SELECT A SUBSET BASED ON INTERSECTION OF 1000 BEST AIC AND BIC MODELS
#=========================================================================

model.list <- lapply(vector.of.models, function(x) lmer(x, data = dataset))
bic_result <- bictab(cand.set = model.list)
aic_result <- aictab(cand.set = model.list)

# pick first 2000 best models from AIC and BIC method, and found 1073 models that's common
aic_bic_dt <- merge(as.data.table(aic_result[1:1000,1:3]), as.data.table(bic_result[1:1000,1:3]), by = c('Modnames','K'))
mod_index <- as.integer(gsub(x = aic_bic_dt$Modnames,pattern = 'Mod',replacement = ''))

#======================================
## STEP3: Fit lmer to the 1073 best models
#=====================================

subset.models <- vector.of.models[mod_index]

list.of.fits <- lapply(subset.models, function(x) {
  
  
  formula    <- as.formula(x)
  print(paste0("running model: ", x)) # progress comment

  temp <- copy(dataset)
  
  fold <- 10
  
  ## assign folds
  temp$folds <- sample(fold)
  
  oos_rmse_list <- numeric()
  is_rmse_list <- numeric()
  for(i in 1:fold){ # "i" is the fold to be left out
    
    model <- lmer(formula, data = temp[folds!=i])
    temp[, paste0("resid_", i)] <- predict(model, newdata=temp,allow.new.levels = T) - temp$y

    # get oos rmse for i th fold, and add to the list
    oos_rmse <- sqrt(mean(temp[folds ==i,get(paste0("resid_", i))]^2, na.rm = T))
    oos_rmse_list[i] <- oos_rmse
    
    # get in sample rmse for i th fold, and add to the list
    is_rmse <- sqrt(mean(temp[folds !=i,get(paste0("resid_", i))]^2, na.rm = T))
    is_rmse_list[i] <- is_rmse
    print(paste0('fold ',i,' has finished'))
  }


  ## then calculate mean rmse from the 10 folds
  oos.rmse <- mean(oos_rmse_list)
  is.rmse <- mean(is_rmse_list)
  
  print('all folds are done, mean RMSE is calculated.')
  
  data.frame(model = x,
             is_rmse = is.rmse,
             oos_rmse  = oos.rmse
             )
})

## Collapse to a data frame, filter top 10 best
result <- as.data.table(do.call(rbind, list.of.fits))
result[order(oos_rmse)][1:10]

