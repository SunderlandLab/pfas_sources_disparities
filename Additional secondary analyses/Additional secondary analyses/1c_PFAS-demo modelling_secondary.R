# ---
# title: 3c_PFAS-demo modelling_secondary
# author: Jahred Liddie
# purpose: regression modeling of relationships between sociodemographic factors (of communities served)
# and PFAS contamination of drinking water; this includes additional county-level race/ethnicity & SES variables
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
                                           0, NA))
                             )

CWS_FN$median.HH.income.thousands <- CWS_FN$median.HH.income/1000

################################################################################
# now, set up iterative modelling
PFAS.vars <- c("PFOA_detect", "PFOS_detect", "PFNA_detect", "PFHxS_detect", "PFBS_detect", "Any_detect")
MCL.vars <- c("above_PFOA_reg", "above_PFOS_reg", "above_any_reg")

f1 <- "~ percHisp + percBlack + poverty + State" # first formula
f2 <- "~ percHisp + percBlack + poverty + percAmind + State" 
f3 <- "~ percHisp + percBlack + poverty + percAPI + State" 
f4 <- "~ percHisp + percBlack + poverty + percAmind + percAPI + State"

f1a <- "~ percHisp + percBlack + poverty + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment" # second formula
f2a <- "~ percHisp + percBlack + poverty + percAmind + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment"
f3a <- "~ percHisp + percBlack + poverty + percAPI + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment"
f4a <- "~ percHisp + percBlack + poverty + percAmind + percAPI + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment"

# this function runs and extracts model output
PFAS_model.f <- function(PFAS.outcome = NULL, model.formula = NULL, dataset = NULL, model.name = NULL) {
  
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
    percent.change = round((exp(coefficient) - 1)*100, 1),
    n.obs = length(PFAS.model$residuals),
    formula = paste(model.formula),
    model = model.name
  )
  
  neat.PFAS.model <- filter(neat.PFAS.model,
                            !grepl("State", dep.var) & dep.var != "(Intercept)"
  )
}

PFAS_models_partial <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f1, dataset = CWS_FN, model.name = "Main model"))
PFAS_models_partial2 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f2, dataset = CWS_FN, model.name = "Model 2 (including % American Indian/Alaskan Native)"))
PFAS_models_partial3 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f3, dataset = CWS_FN, model.name = "Model 3 (including % Asian/Native Hawaiian/other Pacific Islander)"))
PFAS_models_partial4 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f4, dataset = CWS_FN, model.name = "Model 4 (including % American Indian/Alaskan Native & Asian/Native Hawaiian/other Pacific Islander)"))

PFAS_models_fully <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f1a, dataset = CWS_FN, model.name = "Main model"))
PFAS_models_fully2 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f2a, dataset = CWS_FN, model.name = "Model 2 (including % American Indian/Alaskan Native)"))
PFAS_models_fully3 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f3a, dataset = CWS_FN, model.name = "Model 3 (including % Asian/Native Hawaiian/other Pacific Islander)"))
PFAS_models_fully4 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f4a, dataset = CWS_FN, model.name = "Model 4 (including % American Indian/Alaskan Native & Asian/Native Hawaiian/other Pacific Islander)"))

MCL_models_partial <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f1, dataset = CWS_FN, model.name = "Main model"))
MCL_models_partial2 <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f2, dataset = CWS_FN, model.name = "Model 2 (including % American Indian/Alaskan Native)"))
MCL_models_partial3 <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f3, dataset = CWS_FN, model.name = "Model 3 (including % Asian/Native Hawaiian/other Pacific Islander)"))
MCL_models_partial4 <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f4, dataset = CWS_FN, model.name = "Model 4 (including % American Indian/Alaskan Native & Asian/Native Hawaiian/other Pacific Islander)"))

MCL_models_fully <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f1a, dataset = CWS_FN, model.name = "Main model"))
MCL_models_fully2 <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f2a, dataset = CWS_FN, model.name = "Model 2 (including % American Indian/Alaskan Native)"))
MCL_models_fully3 <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f3a, dataset = CWS_FN, model.name = "Model 3 (including % Asian/Native Hawaiian/other Pacific Islander)"))
MCL_models_fully4 <- map_dfr(MCL.vars, ~PFAS_model.f(.x, model.formula = f4a, dataset = CWS_FN, model.name = "Model 4 (including % American Indian/Alaskan Native & Asian/Native Hawaiian/other Pacific Islander)"))

all.PFAS_models_partial <- rbind(PFAS_models_partial, PFAS_models_partial2, PFAS_models_partial3, PFAS_models_partial4)

all.PFAS_models_partial$p.value <- ifelse(round(2*pnorm(q = abs(all.PFAS_models_partial$t.value), lower.tail = F), 3) < 0.001,
                                     "<0.001", round(2*pnorm(q = abs(all.PFAS_models_partial$t.value), lower.tail = F), 3))

all.PFAS_models_fully <- rbind(PFAS_models_fully, PFAS_models_fully2, PFAS_models_fully3, PFAS_models_fully4)

all.PFAS_models_fully$p.value <- ifelse(round(2*pnorm(q = abs(all.PFAS_models_fully$t.value), lower.tail = F), 3) < 0.001,
                                          "<0.001", round(2*pnorm(q = abs(all.PFAS_models_fully$t.value), lower.tail = F), 3))

all.MCL_models_partial <- rbind(MCL_models_partial, MCL_models_partial2, MCL_models_partial3, MCL_models_partial4)

all.MCL_models_partial$p.value <- ifelse(round(2*pnorm(q = abs(all.MCL_models_partial$t.value), lower.tail = F), 3) < 0.001,
                                          "<0.001", round(2*pnorm(q = abs(all.MCL_models_partial$t.value), lower.tail = F), 3))

all.MCL_models_fully <- rbind(MCL_models_fully, MCL_models_fully2, MCL_models_fully3, MCL_models_fully4)

all.MCL_models_fully$p.value <- ifelse(round(2*pnorm(q = abs(all.MCL_models_fully$t.value), lower.tail = F), 3) < 0.001,
                                         "<0.001", round(2*pnorm(q = abs(all.MCL_models_fully$t.value), lower.tail = F), 3))


all.PFAS_models_fully <- all.PFAS_models_fully %>%
  filter(dep.var %in% unique(all.PFAS_models_partial$dep.var))

all.MCL_models_fully <- all.MCL_models_fully %>%
  filter(dep.var %in% unique(all.MCL_models_partial$dep.var))

# tables
reg.table <- all.PFAS_models_partial %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

reg.table2 <- all.PFAS_models_fully %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

reg.table3 <- all.MCL_models_partial %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

reg.table4 <- all.MCL_models_fully %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

# secondary analyses with other SES variables
f5 <- "~ percHisp + percBlack + median.HH.income.thousands + State" 
f6 <- "~ percHisp + percBlack + percunderHSed + State"
f5a <- "~ percHisp + percBlack + median.HH.income.thousands + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment"
f6a <- "~ percHisp + percBlack + percunderHSed + State + any_airport + any_MFTA + 
         WWTP_logflow + any_industry + landfill.LMOP_count + GW_SW + PWS_size_tri + treatment"

PFAS_models_partial5 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f5, dataset = CWS_FN, model.name = "Model 2 (median income)"))
PFAS_models_partial6 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f6, dataset = CWS_FN, model.name = "Model 3 (% with less than HS diploma)"))
PFAS_models_fully5 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f5a, dataset = CWS_FN, model.name = "Model 2 (median income)"))
PFAS_models_fully6 <- map_dfr(PFAS.vars, ~PFAS_model.f(.x, model.formula = f6a, dataset = CWS_FN, model.name = "Model 3 (% with less than HS diploma)"))

all.PFAS_models_partial_SES <- rbind(PFAS_models_partial, PFAS_models_partial5, PFAS_models_partial6)
all.PFAS_models_fully_SES <- rbind(PFAS_models_fully, PFAS_models_fully5, PFAS_models_fully6)
all.PFAS_models_partial_SES$p.value <- ifelse(round(2*pnorm(q = abs(all.PFAS_models_partial_SES$t.value), lower.tail = F), 3) < 0.001,
                                              "<0.001", round(2*pnorm(q = abs(all.PFAS_models_partial_SES$t.value), lower.tail = F), 3))
all.PFAS_models_fully_SES$p.value <- ifelse(round(2*pnorm(q = abs(all.PFAS_models_fully_SES$t.value), lower.tail = F), 3) < 0.001,
                                              "<0.001", round(2*pnorm(q = abs(all.PFAS_models_fully_SES$t.value), lower.tail = F), 3))

all.PFAS_models_fully_SES <- all.PFAS_models_fully_SES %>%
  filter(dep.var %in% unique(all.PFAS_models_partial_SES$dep.var))

# tables
reg.table_SES <- all.PFAS_models_partial_SES %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

reg.table_SES_fully <- all.PFAS_models_fully_SES %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

if (FALSE) {
  write.csv(reg.table, "../../Regressions/PFAS_demographics_secondary1.csv")
  write.csv(reg.table2, "../../Regressions/PFAS_demographics_secondary2.csv")
  write.csv(reg.table3, "../../Regressions/PFAS_demographics_secondary3.csv")
  write.csv(reg.table4, "../../Regressions/PFAS_demographics_secondary4.csv")
  
  write.csv(reg.table_SES, "../../Regressions/PFAS_SES_secondary.csv")
  write.csv(reg.table_SES_fully, "../../Regressions/PFAS_SES_fulladjust_secondary.csv")
  
}

