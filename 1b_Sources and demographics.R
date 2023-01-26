# ---
# title: 3b_Sources and demographics
# author: Jahred Liddie
# purpose: regression modeling of relationships between sociodemographic factors (of communities served)
  # and contamination sources
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

CWS_FN <- CWS_FN %>%
  mutate(across(c(any_airport, any_MFTA, any_industry), ~ifelse(.x == "One or more", 1, 0)))

contam.vars <- c("any_MFTA", "any_airport", "any_industry", "WWTP_flow.median", "landfill_median")

# first formula: need to use this fixed effect notation for conleyreg
f1 <- "~ percHisp + percBlack + poverty + State"

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

state.subset <- subset(CWS_FN, State != "CA" & State != "PA" & State != "MD" & State != "UT")

select.states_models <- map_dfr(contam.vars, ~sources_model.f(.x, model.formula = f1, dataset = state.subset))

all.sources_models$p.value <- ifelse(round(2*pnorm(q = abs(all.sources_models$t.value), lower.tail = F), 3) < 0.001,
                                     "<0.001", round(2*pnorm(q = abs(all.sources_models$t.value), lower.tail = F), 3))

all.sources_models <- all.sources_models %>%
  group_by(outcome) %>%
  arrange(desc(coefficient)) %>%
  mutate(coefficient.order = row_number(n.obs)) %>%
  ungroup()

# for regression table
reg.table <- all.sources_models %>%
  mutate(table.coef = paste(round(percent.change, 2), 
                            ifelse(p.value == "<0.001" | as.numeric(p.value) < 0.01, "***",
                                   ifelse(as.numeric(p.value) < 0.05, "**",
                                          ifelse(as.numeric(p.value) < 0.1, "*", ""))),
                            " (", LCI, ", " , UCI, ")", sep = "")
         ) %>%
  pivot_wider(id_cols = "dep.var",
              names_from = "outcome",
              values_from = "table.coef")

if (FALSE) {
  write.csv(reg.table, "../Regressions/sources_demographics.csv")
}
