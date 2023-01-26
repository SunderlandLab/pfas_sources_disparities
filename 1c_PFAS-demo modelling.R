# ---
# title: 3c_PFAS-demo modelling
# author: Jahred Liddie
# purpose: regression modeling of relationships between sociodemographic factors (of communities served)
# and PFAS contamination of drinking water
# date: 7/30/2022
# ---

library(tidyverse)
library(sandwich)
library(lmtest)
library(splines)
library(mgcv)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# note: these datasets are loaded (and already processed for analysis) from the statewide sampling repository
CWS_FN <- read.csv("/Users/jahredl/Dropbox/State PFAS data/Merged PFAS data/final_dat_FN.csv")
CWS_RW <- read.csv("/Users/jahredl/Dropbox/State PFAS data/Merged PFAS data/final_dat_RW.csv")

################################################################################
## DEFINE FUNCTIONS ##
# create CRSE function for logistic models
logistic.f <- function(dataframe, formula, cluster) {
    
    model <- glm(formula, data = dataframe, family = "binomial")
  
    # adjustment for clustering
    crse <- coeftest(model, vcov. = vcovCL(model, cluster = cluster))
    
    varcov <- vcovCL(model, cluster = cluster) # save var/covar matrix
    
    return(list(model, crse, varcov))
}

################################################################################
# first change this to a factor and re-level to NJ
CWS_FN$State <- as.factor(CWS_FN$State)
CWS_FN$State <- relevel(CWS_FN$State, ref = "NJ")

CWS_RW$State <- as.factor(CWS_RW$State)
CWS_RW$State <- relevel(CWS_RW$State, ref = "NJ")
# logarithm of total existing effluent; if flow == 0, replace with 0 for log flow
CWS_FN$WWTP_logflow <- ifelse(CWS_FN$WWTP_existingtotalflow_Mgal.d == 0 , 0,
                           log(CWS_FN$WWTP_existingtotalflow_Mgal.d))
CWS_RW$WWTP_logflow <- ifelse(CWS_RW$WWTP_existingtotalflow_Mgal.d == 0 , 0,
                              log(CWS_RW$WWTP_existingtotalflow_Mgal.d))

################################################################################
# model exploration prior to iterative process
################################################################################
# models without covariates or interactions
n1 <- logistic.f(dataframe = CWS_FN, 
                      formula = "PFOA_detect ~ percHisp + percBlack + poverty + State",
                      cluster = CWS_FN$state_county)

# investigating at potential non-linearity
n1.gam <- gam(PFOA_detect ~ s(percHisp, bs = "cr", fx = F) + 
                s(percBlack, bs = "cr", fx = F) + 
                s(poverty, bs = "cr", fx = F) + State,
              data = CWS_FN)

# adjustment for clustering
crse <- coeftest(n1.gam, vcov. = vcovCL(n1.gam, cluster = CWS_FN$state_county))

# generally, after adjustment for county-level clustering, no evidence of non-linearity (see DOFs in crse)

# including additonal 
o1 <- logistic.f(dataframe = CWS_FN, 
                 formula = "PFOA_detect ~ percHisp + percBlack + poverty + State + 
                                          any_airport + any_MFTA +  WWTP_logflow + any_industry +
                                          GW_SW + PWS_size_tri + treatment",
                 cluster = CWS_FN$state_county)

# add any reg indicator
CWS_FN$above_any_reg <- with(data = CWS_FN,
                             ifelse((above_PFOA_reg > 0 |
                                       above_PFOS_reg > 0 |
                                       above_PFNA_reg > 0 | 
                                       above_PFHxS_reg > 0 |
                                       above_PFBS_reg > 0) &
                                      (!is.na(above_PFOA_reg) & !is.na(above_PFOS_reg) &
                                         !is.na(above_PFNA_reg) & !is.na(above_PFHxS_reg) &
                                         !is.na(above_PFBS_reg)), 1,
                                    ifelse(!is.na(above_PFOA_reg) & !is.na(above_PFOS_reg) &
                                             !is.na(above_PFNA_reg) & !is.na(above_PFHxS_reg) &
                                             !is.na(above_PFBS_reg),
                                           0, NA)))

MCL_all <- logistic.f(CWS_FN,
                      formula = "above_any_reg ~ percHisp + percBlack + poverty + State",
                      cluster = CWS_FN$state_county)

MCL_all_sources <- logistic.f(dataframe = CWS_FN,
                              formula = "above_any_reg ~ percHisp + percBlack + poverty + GW_SW + PWS_size_tri +
                               any_airport + any_MFTA + landfill.LMOP_count + treatment + any_industry + State",
                              cluster = CWS_FN$state_county)

################################################################################
# now, set up iterative modelling
PFAS.vars <- c("PFOA_detect", "PFOS_detect", "PFNA_detect", "PFHxS_detect", "PFBS_detect", "Any_detect")
MCL.vars <- c("above_PFOA_reg", "above_PFOS_reg", "above_any_reg")
f1 <- "~ percHisp + percBlack + poverty + State" # first formula
f2 <- "~ percHisp + percBlack + poverty + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment" # second formula

# this function runs and extracts model output
PFAS_model.f <- function(PFAS.outcome = NULL, model.formula = NULL, dataset = NULL) {
  
  PFAS.model <- glm(paste(PFAS.outcome, model.formula),
                       family = "binomial",
                       data = dataset)
  
  # cluster robust standard errors (county as the cluster)
  crse <- coeftest(PFAS.model, vcov. = vcovCL(PFAS.model, cluster = dataset$state_county))
  
  SE <- tibble(
    dep.var = names(PFAS.model$coefficients),
    SE = crse[, "Std. Error"]
  )
  
  # summary table
  neat.PFAS.model <- tibble(
    outcome = rep(PFAS.outcome, length(coef(PFAS.model))),
    # outcome.name = str_replace(outcome, pattern = "log", replacement = ""),
    dep.var = names(coef(PFAS.model)),
    coefficient = coef(PFAS.model),
    OR = exp(coefficient),
    coef.se = SE$SE[mapply(grepl, dep.var, SE$dep.var)],
    UCI.OR = exp(coefficient + 1.96*coef.se),
    LCI.OR = exp(coefficient - 1.96*coef.se),
    UCI = round((exp(coefficient + 1.96*coef.se) - 1)*100, 1),
    LCI = round((exp(coefficient - 1.96*coef.se) - 1)*100, 1),
    t.value = coefficient/coef.se,
    p.value = 2*pnorm(q = abs(t.value), lower.tail = FALSE),
    # using z approximation for p-values and CIs (large sample size)
    percent.change = paste(round((exp(coefficient) - 1)*100, 1), # formatted for final tables
                           ifelse(p.value < 0.01, "***", 
                                  ifelse(p.value < 0.05, "**",
                                         ifelse(p.value < 0.1, "*", ""))),
                           sep = ""),
    n.obs = length(PFAS.model$residuals),
    formula = paste(model.formula)
  )
  
  neat.PFAS.model <- filter(neat.PFAS.model,
                            !grepl("State", dep.var) & dep.var != "(Intercept)"
  )
}

PFAS_models_partial <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f1, dataset = CWS_FN))

PFAS_models_fully <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f2, dataset = CWS_FN))

MCL_models_partial <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f1, dataset = CWS_FN))

MCL_models_fully <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f2, dataset = CWS_FN))

################################################################################
#### Sensitivity analysis:
# this uses a *slightly* different dataset in which samples that were not clearly denoted as
# raw or finished water were assumed to be raw water (which impacts data processing)
################################################################################
# only for comparison with other dataset
o1.RW <- logistic.f(dataframe = CWS_RW, 
                    formula = "PFOA_detect ~ percHisp + percBlack + poverty + State + 
                                          any_airport + any_MFTA +  WWTP_logflow + any_industry +
                                          GW_SW + PWS_size_tri + treatment",
                    cluster = CWS_RW$state_county)

o2.RW <- logistic.f(dataframe = CWS_RW, 
                    formula = "PFOS_detect ~ percHisp + percBlack + poverty + State + 
                                          any_airport + any_MFTA +  WWTP_logflow + any_industry +
                                          GW_SW + PWS_size_tri + treatment",
                    cluster = CWS_RW$state_county)

o6.RW <- logistic.f(dataframe = CWS_RW, 
                    formula = "Any_detect ~ percHisp + percBlack + poverty + State +
                                         any_airport + any_MFTA +  WWTP_logflow + any_industry +
                                         GW_SW + PWS_size_tri + treatment",
                    cluster = CWS_RW$state_county)

# results very similar compared to prior results 
################################################################################
rm(list = c(ls(pattern = "RW")))
