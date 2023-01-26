# ---
# title: 3d_Secondary analysis modelling
# purpose: secondary analyses (and one additional sensitivity analysis)
# date: 7/30/2022
# ---

library(tidyverse)
library(sandwich)
library(lmtest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# note: this dataset is loaded (and already processed for analysis) from the statewide sampling repository
CWS_FN <- read.csv("/Users/jahredl/Dropbox/State PFAS data/Merged PFAS data/final_dat_FN.csv")

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

# logarithmic flow; if flow == 0, replace with 0 for log flow
CWS_FN$WWTP_logflow <- ifelse(CWS_FN$WWTP_existingtotalflow_Mgal.d == 0 , 0,
                           log(CWS_FN$WWTP_existingtotalflow_Mgal.d))

CWS_FN$WWTP_logflow.p.area <- ifelse(CWS_FN$WWTP_existingtotalflow_Mgal.d.sqmi == 0 , 0,
                                  log(CWS_FN$WWTP_existingtotalflow_Mgal.d.sqmi))

CWS_FN$landfill_log.count <- ifelse(CWS_FN$landfill.LMOP_count == 0 , 0,
                                     log(CWS_FN$landfill.LMOP_count))

################################################################################
# secondary analysis: flow vs flow per area for WWTPs - any differences in estimates of interest?
  # also adding count of landfills to main fully adjusted model

# also do these after filtering out states with specific sampling strategies
CWS_excl <- CWS_FN %>% filter(State != "CA" & State != "PA" & State != "MD" & State != "NY")

# now, set up iterative modelling
PFAS.vars <- c("PFOA_detect", "PFOS_detect", "Any_detect")
f1 <- "~ percHisp + percBlack + poverty + State" # main partially adjusted formula
f2 <- "~ percHisp + percBlack + poverty + State + any_airport + 
         any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment" # main fully adjusted formula
f3 <- "~ percHisp + percBlack + poverty + State + any_airport + 
         any_MFTA + WWTP_logflow.p.area + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment" 
        # secondary analysis formula (flow per area of HUC)

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

PFAS_models_fully_sec <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f2, dataset = CWS_FN))

PFAS_models_sec1 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f3, dataset = CWS_FN))

################################################################################
# stratifying by urban/rural status
CWS_FN$urban <- ifelse(CWS_FN$percurban >= 50, 1, 0)

PFAS_models_rural_partial <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f1, 
                                                              dataset = subset(CWS_FN, urban == 0)))
PFAS_models_urban_partial <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f1, 
                                                              dataset = subset(CWS_FN, urban == 1)))

PFAS_models_rural_fully <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f2, 
                                                            dataset = subset(CWS_FN, urban == 0)))
PFAS_models_urban_fully <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f2, 
                                                            dataset = subset(CWS_FN, urban == 1)))

secondary.results <- rbind(PFAS_models_urban_partial %>% mutate(type = "Stratified: urban (partial)"),
                           PFAS_models_urban_fully %>% mutate(type = "Stratified: urban (fully)"),
                           PFAS_models_rural_partial %>% mutate(type = "Stratified: rural (partial)"),
                           PFAS_models_rural_fully %>% mutate(type = "Stratified: rural (fully)"))

write.csv(secondary.results, "../Merged PFAS data/Results/secondary_stratified_results.csv")

################################################################################
# sensitivity using more system categories
f_size <- "~ percHisp + percBlack + poverty + State + 
             any_airport + any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + 
             GW_SW + PWS_size + treatment"

PFAS_models_moresize <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f_size, dataset = CWS_FN))

