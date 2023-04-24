# ---
# title: Secondary analysis modelling
# author: Jahred Liddie
# purpose: secondary analyses and additional sensitivity analyses
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
    
    n.obs <- length(model$residuals)
    
    return(list(model, crse, varcov, n.obs))
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

if (FALSE) {
write.csv(secondary.results, "../Merged PFAS data/Results/secondary_stratified_results.csv")
}

################################################################################
# sensitivity using more system categories
f_size <- "~ percHisp + percBlack + poverty + State + 
             any_airport + any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + 
             GW_SW + PWS_size + treatment"

PFAS_models_moresize <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f_size, dataset = CWS_FN))

################################################################################
# sensitivity / bounds analysis substituting missing data for Any_detect
# this uses a simple substitution where: any system reporting detection of ≥1 PFAS,
# but missing ≥1 PFAS gets Any_detect == 1. If there are no detections of any compounds but
# missing ≥1 PFAS, assuming that the system had Any_detect == 0.
CWS_FN <- CWS_FN %>%
  mutate(Any_detect2 = case_when(Any_detect == 0 ~ 0,
                                 Any_detect == 1 ~ 1,
                                 is.na(Any_detect) & (PFOA_detect == 1 |
                                                      PFNA_detect == 1 |
                                                      PFOS_detect == 1 |
                                                      PFBS_detect == 1 |
                                                      PFHxS_detect == 1) ~ 1,
                                 TRUE ~ 0)
         )

n1_allAny_detect <- logistic.f(formula = "Any_detect2 ~ percHisp + percBlack + poverty + State",
                               dataframe = CWS_FN, cluster = CWS_FN$state_county)

o1_allAny_detect <- logistic.f(formula = "Any_detect2 ~ percHisp + percBlack + poverty + State +
                                         any_airport + any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + 
                                         GW_SW + PWS_size + treatment", dataframe = CWS_FN,
                               cluster = CWS_FN$state_county)

# main models
n1_Any_detect <- logistic.f(formula = "Any_detect ~ percHisp + percBlack + poverty + State",
                               dataframe = CWS_FN, cluster = CWS_FN$state_county)

o1_Any_detect <- logistic.f(formula = "Any_detect ~ percHisp + percBlack + poverty + State +
                                         any_airport + any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + 
                                         GW_SW + PWS_size + treatment", dataframe = CWS_FN,
                            cluster = CWS_FN$state_county)

missing.results <- tibble(
  outcome = c("≥ 1 PFAS (of 5 total) [all cases, partial]",
              "≥ 1 PFAS (of 5 total) [all cases, full]",
              "≥ 1 PFAS (of 5 total) [complete cases, partial]",
              "≥ 1 PFAS (of 5 total) [complete cases, full]"),
  percHisp.odds = c(n1_allAny_detect[[2]]["percHisp","Estimate"], 
                    o1_allAny_detect[[2]]["percHisp","Estimate"], 
                    n1_Any_detect[[2]]["percHisp","Estimate"], 
                    o1_Any_detect[[2]]["percHisp","Estimate"]),
  percHisp.sd = c(n1_allAny_detect[[2]]["percHisp","Std. Error"], 
                  o1_allAny_detect[[2]]["percHisp","Std. Error"], 
                  n1_Any_detect[[2]]["percHisp","Std. Error"], 
                  o1_Any_detect[[2]]["percHisp","Std. Error"]),
  percBlack.odds = c(n1_allAny_detect[[2]]["percBlack","Estimate"], 
                     o1_allAny_detect[[2]]["percBlack","Estimate"], 
                     n1_Any_detect[[2]]["percBlack","Estimate"], 
                     o1_Any_detect[[2]]["percBlack","Estimate"]),
  percBlack.sd = c(n1_allAny_detect[[2]]["percBlack","Std. Error"], 
                   o1_allAny_detect[[2]]["percBlack","Std. Error"], 
                   n1_Any_detect[[2]]["percBlack","Std. Error"], 
                   o1_Any_detect[[2]]["percBlack","Std. Error"]),
  poverty.odds = c(n1_allAny_detect[[2]]["poverty","Estimate"], 
                   o1_allAny_detect[[2]]["poverty","Estimate"], 
                   n1_Any_detect[[2]]["poverty","Estimate"], 
                   o1_Any_detect[[2]]["poverty","Estimate"]),
  poverty.sd = c(n1_allAny_detect[[2]]["poverty","Std. Error"], 
                 o1_allAny_detect[[2]]["poverty","Std. Error"], 
                 n1_Any_detect[[2]]["poverty","Std. Error"], 
                 o1_Any_detect[[2]]["poverty","Std. Error"]),
  percHisp.p.value = c(n1_allAny_detect[[2]]["percHisp","Pr(>|z|)"], 
                       o1_allAny_detect[[2]]["percHisp","Pr(>|z|)"], 
                       n1_Any_detect[[2]]["percHisp","Pr(>|z|)"], 
                       o1_Any_detect[[2]]["percHisp","Pr(>|z|)"]),
  percBlack.p.value = c(n1_allAny_detect[[2]]["percBlack","Pr(>|z|)"], 
                        o1_allAny_detect[[2]]["percBlack","Pr(>|z|)"], 
                        n1_Any_detect[[2]]["percBlack","Pr(>|z|)"], 
                        o1_Any_detect[[2]]["percBlack","Pr(>|z|)"]),
  poverty.p.value = c(n1_allAny_detect[[2]]["poverty","Pr(>|z|)"], 
                      o1_allAny_detect[[2]]["poverty","Pr(>|z|)"], 
                      n1_Any_detect[[2]]["poverty","Pr(>|z|)"], 
                      o1_Any_detect[[2]]["poverty","Pr(>|z|)"]),
  percHisp.full.estimate = paste(round((exp(percHisp.odds)-1)*100, 1), 
                                 " [95% CI: ", round((exp(percHisp.odds - 1.96*percHisp.sd)-1)*100, 1), ", ",
                                 round((exp(percHisp.odds + 1.96*percHisp.sd)-1)*100, 1), "]", sep = ""),
  percBlack.full.estimate = paste(round((exp(percBlack.odds)-1)*100, 1), 
                                 " [95% CI: ", round((exp(percBlack.odds - 1.96*percBlack.sd)-1)*100, 1), ", ",
                                 round((exp(percBlack.odds + 1.96*percBlack.sd)-1)*100, 1), "]", sep = ""),
  poverty.full.estimate = paste(round((exp(poverty.odds)-1)*100, 1), 
                                 " [95% CI: ", round((exp(poverty.odds - 1.96*poverty.sd)-1)*100, 1), ", ",
                                 round((exp(poverty.odds + 1.96*poverty.sd)-1)*100, 1), "]", sep = ""),
  n.obs = c(n1_allAny_detect[[4]],
            o1_allAny_detect[[4]],
            n1_Any_detect[[4]],
            o1_Any_detect[[4]])
)

missing.results <- missing.results %>%
  mutate(percHisp.full.estimate = ifelse(percHisp.p.value < 0.01, paste(percHisp.full.estimate, "***", sep = ""),
                                         ifelse(percHisp.p.value < 0.05, paste(percHisp.full.estimate, "**", sep = ""),
                                                ifelse(percHisp.p.value < 0.1, paste(percHisp.full.estimate, "*", sep = ""), percHisp.full.estimate))),
         percBlack.full.estimate = ifelse(percBlack.p.value < 0.01, paste(percBlack.full.estimate, "***", sep = ""),
                                          ifelse(percBlack.p.value < 0.05, paste(percBlack.full.estimate, "**", sep = ""),
                                                 ifelse(percBlack.p.value < 0.1, paste(percBlack.full.estimate, "*", sep = ""), percBlack.full.estimate))),
         poverty.full.estimate = ifelse(poverty.p.value < 0.01, paste(poverty.full.estimate, "***", sep = ""),
                                        ifelse(poverty.p.value < 0.05, paste(poverty.full.estimate, "**", sep = ""),
                                               ifelse(poverty.p.value < 0.1, paste(poverty.full.estimate, "*", sep = ""), poverty.full.estimate)))
         )

missing.results.table <- t(missing.results %>% select(outcome, ends_with("estimate"), n.obs))

if (FALSE) {
  write.csv(missing.results.table, "../Regressions/missing_regressions.csv")
}

# stratified models
n1_allAny_detect_urban <- logistic.f(formula = "Any_detect2 ~ percHisp + percBlack + poverty + State",
                               dataframe = CWS_FN %>% filter(urban == 1), cluster = CWS_FN$state_county[CWS_FN$urban==1])

o1_allAny_detect_urban <- logistic.f(formula = "Any_detect2 ~ percHisp + percBlack + poverty + State +
                                         any_airport + any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + 
                                         GW_SW + PWS_size + treatment", dataframe = CWS_FN %>% filter(urban == 1),
                               cluster = CWS_FN$state_county[CWS_FN$urban==1])

n1_allAny_detect_rural <- logistic.f(formula = "Any_detect2 ~ percHisp + percBlack + poverty + State",
                                     dataframe = CWS_FN %>% filter(urban == 0), cluster = CWS_FN$state_county[CWS_FN$urban==0])

o1_allAny_detect_rural <- logistic.f(formula = "Any_detect2 ~ percHisp + percBlack + poverty + State +
                                         any_airport + any_MFTA + WWTP_logflow + any_industry + landfill.LMOP_count + 
                                         GW_SW + PWS_size + treatment", dataframe = CWS_FN %>% filter(urban == 0),
                                     cluster = CWS_FN$state_county[CWS_FN$urban==0])
missing.results.stratified <- tibble(
  outcome = c("≥ 1 PFAS (of 5 total) [complete cases, partial] - urban",
              "≥ 1 PFAS (of 5 total) [complete cases, full] - urban",
              "≥ 1 PFAS (of 5 total) [complete cases, partial] - rural",
              "≥ 1 PFAS (of 5 total) [complete cases, full] - rural"),
  percHisp.odds = c(n1_allAny_detect_urban[[2]]["percHisp","Estimate"], 
                    o1_allAny_detect_urban[[2]]["percHisp","Estimate"], 
                    n1_allAny_detect_rural[[2]]["percHisp","Estimate"], 
                    o1_allAny_detect_rural[[2]]["percHisp","Estimate"]),
  percHisp.sd = c(n1_allAny_detect_urban[[2]]["percHisp","Std. Error"], 
                  o1_allAny_detect_urban[[2]]["percHisp","Std. Error"], 
                  n1_allAny_detect_rural[[2]]["percHisp","Std. Error"], 
                  o1_allAny_detect_rural[[2]]["percHisp","Std. Error"]),
  percBlack.odds = c(n1_allAny_detect_urban[[2]]["percBlack","Estimate"], 
                     o1_allAny_detect_urban[[2]]["percBlack","Estimate"], 
                     n1_allAny_detect_rural[[2]]["percBlack","Estimate"], 
                     o1_allAny_detect_rural[[2]]["percBlack","Estimate"]),
  percBlack.sd = c(n1_allAny_detect_urban[[2]]["percBlack","Std. Error"], 
                   o1_allAny_detect_urban[[2]]["percBlack","Std. Error"], 
                   n1_allAny_detect_rural[[2]]["percBlack","Std. Error"], 
                   o1_allAny_detect_rural[[2]]["percBlack","Std. Error"]),
  poverty.odds = c(n1_allAny_detect_urban[[2]]["poverty","Estimate"], 
                   o1_allAny_detect_urban[[2]]["poverty","Estimate"], 
                   n1_allAny_detect_rural[[2]]["poverty","Estimate"], 
                   o1_allAny_detect_rural[[2]]["poverty","Estimate"]),
  poverty.sd = c(n1_allAny_detect_urban[[2]]["poverty","Std. Error"], 
                 o1_allAny_detect_urban[[2]]["poverty","Std. Error"], 
                 n1_allAny_detect_rural[[2]]["poverty","Std. Error"], 
                 o1_allAny_detect_rural[[2]]["poverty","Std. Error"]),
  percHisp.p.value = c(n1_allAny_detect_urban[[2]]["percHisp","Pr(>|z|)"], 
                       o1_allAny_detect_urban[[2]]["percHisp","Pr(>|z|)"], 
                       n1_allAny_detect_rural[[2]]["percHisp","Pr(>|z|)"], 
                       o1_allAny_detect_rural[[2]]["percHisp","Pr(>|z|)"]),
  percBlack.p.value = c(n1_allAny_detect_urban[[2]]["percBlack","Pr(>|z|)"], 
                        o1_allAny_detect_urban[[2]]["percBlack","Pr(>|z|)"], 
                        n1_allAny_detect_rural[[2]]["percBlack","Pr(>|z|)"], 
                        o1_allAny_detect_rural[[2]]["percBlack","Pr(>|z|)"]),
  poverty.p.value = c(n1_allAny_detect_urban[[2]]["poverty","Pr(>|z|)"], 
                      o1_allAny_detect_urban[[2]]["poverty","Pr(>|z|)"], 
                      n1_allAny_detect_rural[[2]]["poverty","Pr(>|z|)"], 
                      o1_allAny_detect_rural[[2]]["poverty","Pr(>|z|)"]),
  percHisp.full.estimate = paste(round((exp(percHisp.odds)-1)*100, 1), 
                                 " [95% CI: ", round((exp(percHisp.odds - 1.96*percHisp.sd)-1)*100, 1), ", ",
                                 round((exp(percHisp.odds + 1.96*percHisp.sd)-1)*100, 1), "]", sep = ""),
  percBlack.full.estimate = paste(round((exp(percBlack.odds)-1)*100, 1), 
                                  " [95% CI: ", round((exp(percBlack.odds - 1.96*percBlack.sd)-1)*100, 1), ", ",
                                  round((exp(percBlack.odds + 1.96*percBlack.sd)-1)*100, 1), "]", sep = ""),
  poverty.full.estimate = paste(round((exp(poverty.odds)-1)*100, 1), 
                                " [95% CI: ", round((exp(poverty.odds - 1.96*poverty.sd)-1)*100, 1), ", ",
                                round((exp(poverty.odds + 1.96*poverty.sd)-1)*100, 1), "]", sep = ""),
  n.obs = c(n1_allAny_detect_urban[[4]],
            o1_allAny_detect_urban[[4]],
            n1_allAny_detect_rural[[4]],
            o1_allAny_detect_rural[[4]])
)

missing.results.stratified <- missing.results.stratified %>%
  mutate(percHisp.full.estimate = ifelse(percHisp.p.value < 0.01, paste(percHisp.full.estimate, "***", sep = ""),
                                         ifelse(percHisp.p.value < 0.05, paste(percHisp.full.estimate, "**", sep = ""),
                                                ifelse(percHisp.p.value < 0.1, paste(percHisp.full.estimate, "*", sep = ""), percHisp.full.estimate))),
         percBlack.full.estimate = ifelse(percBlack.p.value < 0.01, paste(percBlack.full.estimate, "***", sep = ""),
                                          ifelse(percBlack.p.value < 0.05, paste(percBlack.full.estimate, "**", sep = ""),
                                                 ifelse(percBlack.p.value < 0.1, paste(percBlack.full.estimate, "*", sep = ""), percBlack.full.estimate))),
         poverty.full.estimate = ifelse(poverty.p.value < 0.01, paste(poverty.full.estimate, "***", sep = ""),
                                        ifelse(poverty.p.value < 0.05, paste(poverty.full.estimate, "**", sep = ""),
                                               ifelse(poverty.p.value < 0.1, paste(poverty.full.estimate, "*", sep = ""), poverty.full.estimate)))
  )

missing.results.stratified.table <- t(missing.results.stratified %>% select(outcome, ends_with("estimate"), n.obs))

if (FALSE) {
  write.csv(missing.results.stratified.table, "../Regressions/missing_regressions_stratified.csv")
}
