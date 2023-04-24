# ---
# title: Sources and demographics_secondary
# author: Jahred Liddie
# purpose: regression modeling of relationships between sociodemographic factors (of communities served)
# and contamination sources; this includes additional county-level race/ethnicity and SES variables
# date: 7/30/2022
# ---

library(tidyverse)
library(sandwich)
library(lmtest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# note: this dataset is loaded (and already processed for analysis) from the statewide sampling repository
CWS_FN <- read.csv("/Users/jahredl/Dropbox/State PFAS data/Merged PFAS data/final_dat_FN.csv")

CWS_FN$WWTP_flow.median <- ifelse(CWS_FN$WWTP_existingtotalflow_Mgal.d <=
                                    median(CWS_FN$WWTP_existingtotalflow_Mgal.d),
                                  0, 1)

CWS_FN$landfill_median <- ifelse(CWS_FN$landfill.LMOP_count <=
                                   median(CWS_FN$landfill.LMOP_count),
                                 0, 1)

CWS_FN$median.HH.income.thousands <- CWS_FN$median.HH.income/1000

CWS_FN <- CWS_FN %>%
  mutate(across(c(any_airport, any_MFTA, any_industry), ~ifelse(.x == "One or more", 1, 0)))

contam.vars <- c("any_MFTA", "any_airport", "any_industry", "WWTP_flow.median", "landfill_median")

# formulas
f1 <- "~ percHisp + percBlack + poverty + State"
f2 <- "~ percHisp + percBlack + poverty + percAmind + State"
f3 <- "~ percHisp + percBlack + poverty + percAPI + State"
f4 <- "~ percHisp + percBlack + poverty + percAPI + percAmind + State"

f5 <- "~ percHisp + percBlack + median.HH.income.thousands + State"
f6 <- "~ percHisp + percBlack + percunderHSed + State"

# note that variance of percAPI and percAmind are much lower than other race/ethnicity vars
# 11 - 36% and 0.8 - 2.7%, respectively. so, may need caution in using these variables, particularly when stratifying
var(CWS_FN$percAmind) / var(CWS_FN$percBlack) *100
var(CWS_FN$percAmind) / var(CWS_FN$percHisp) *100

var(CWS_FN$percAPI) / var(CWS_FN$percBlack) *100
var(CWS_FN$percAPI) / var(CWS_FN$percHisp) *100

# check prior to modeling: how many sum to 100% when including all of these 4?
CWS_FN <- CWS_FN %>% 
  rowwise() %>%
  mutate(flag_sum = sum(percAmind + percBlack + percAPI + percHisp))

# this function runs and extracts model output
sources_model.f <- function(contam.outcome = NULL, model.formula = NULL, dataset = NULL) {
  
  sources.model <- glm(paste(contam.outcome, model.formula),
                       family = "binomial",
                       data = dataset)
  
  crse <- coeftest(sources.model, vcov. = vcovCL(sources.model, cluster = dataset$state_county))
  
  SE <- tibble(
    dep.var = names(sources.model$coefficients),
    SE = crse[, "Std. Error"]
  )
  
  # summary table
  neat.sources.model <- tibble(
    outcome = rep(contam.outcome, length(coef(sources.model))),
    # outcome.name = str_replace(outcome, pattern = "log", replacement = ""),
    dep.var = names(coef(sources.model)),
    coefficient = coef(sources.model),
    OR = exp(coefficient),
    percent.change = round((exp(coefficient) - 1)*100, 1),
    coef.se = SE$SE[mapply(grepl, dep.var, SE$dep.var)],
    UCI.OR = exp(coefficient + 1.96*coef.se),
    LCI.OR = exp(coefficient - 1.96*coef.se),
    UCI = round((exp(coefficient + 1.96*coef.se) - 1)*100, 1),
    LCI = round((exp(coefficient - 1.96*coef.se) - 1)*100, 1),
    t.value = coefficient/coef.se,
    # using z approximation
    n.obs = length(sources.model$residuals),
    formula = paste(model.formula)
  )
  
  neat.sources.model <- filter(neat.sources.model,
                               !grepl("State", dep.var) & dep.var != "(Intercept)"
  )
}

all.sources_models <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f1, dataset = CWS_FN))
all.sources_models2 <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f2, dataset = CWS_FN))
all.sources_models3 <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f3, dataset = CWS_FN))
all.sources_models4 <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f4, dataset = CWS_FN))

all.sources_models$model <- "Main model"
all.sources_models2$model <- "Model 2 (including % American Indian/Alaskan Native)"
all.sources_models3$model <- "Model 3 (including % Asian/Native Hawaiian/other Pacific Islander)"
all.sources_models4$model <- "Model 4 (including % American Indian/Alaskan Native & Asian/Native Hawaiian/other Pacific Islander)"

all.sources_models_demo <- rbind(all.sources_models, all.sources_models2, all.sources_models3, all.sources_models4)

all.sources_models_demo$p.value <- ifelse(round(2*pnorm(q = abs(all.sources_models_demo$t.value), lower.tail = F), 3) < 0.001,
                                     "<0.001", round(2*pnorm(q = abs(all.sources_models_demo$t.value), lower.tail = F), 3))

# for regression table
reg.table_demo <- all.sources_models_demo %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
  ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = c("outcome","model"),
              values_from = "table.coef")

all.sources_models5 <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f5, dataset = CWS_FN))
all.sources_models6 <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f6, dataset = CWS_FN))

all.sources_models5$model <- "Model 2 (median income)"
all.sources_models6$model <- "Model 4 (% with less than HS diploma)"

all.sources_models_SES <- rbind(all.sources_models, all.sources_models5, all.sources_models6)

all.sources_models_SES$p.value <- ifelse(round(2*pnorm(q = abs(all.sources_models_SES$t.value), lower.tail = F), 3) < 0.001,
                                     "<0.001", round(2*pnorm(q = abs(all.sources_models_SES$t.value), lower.tail = F), 3))

reg.table_SES <- all.sources_models_SES %>%
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
  
  write.csv(reg.table_demo, "../../Regressions/sources_demographics_secondary.csv")
  write.csv(reg.table_SES, "../../Regressions/sources_SES_secondary.csv")
  
}