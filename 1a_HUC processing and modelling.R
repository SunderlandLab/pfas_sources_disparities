# ---
# title: HUC processing and modelling
# author: Jahred Liddie
# purpose: data processing; regression modeling of relationships between contamination sources and 
  # PFAS concentrations in 8-digit HUCs
# date: 7/30/2022
# ---

library(tidyverse)
library(spdep)
library(spatialreg)
library(sf)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

HUCs <- st_read("../PFAS point source data/huc250k_shp/huc250k_shp/huc250k.shp")

# note: these two datasets are loaded (and already processed for analysis) from the statewide sampling repository
CWS <- read.csv("/Users/jahredl/Dropbox/State PFAS data/Merged PFAS data/final_dat_FN.csv")
all.measurements <- read.csv("/Users/jahredl/Dropbox/State PFAS data/Merged PFAS data/dat.csv")

all.measurements <- all.measurements %>% filter(PWSID %in% CWS$PWSID)

HUCs <- HUCs %>% select(-HUC250K_, -HUC250K_ID)

## point source data ##
# military sites: this is the EWG list
military <- geojsonsf::geojson_sf("../PFAS point source data/Military_2020MARCH04.geojson")
military <- st_transform(military, crs = st_crs(HUCs)) # convert to NAD83

# suspected <- geojsonio::geojson_read("PFAS point source data/suspected_sites_2020MARCH23.geojson", what = "sp")
# suspected <- st_as_sf(suspected)
# suspected <- st_transform(suspected, crs = "NAD83") # convert to NAD83

# Part 139 certified airports
airports <- readxl::read_excel("../PFAS point source data/Part 139_cert_airports.xlsx")
airports <- as.data.frame(airports)

# EPA PFOS/PFOA stewardship list
epa <- readxl::read_excel("../PFAS point source data/EPA 2010.2015 PFOA Stewardship Program sites.xlsx")

# CWNS wastewater treatment plants
WWTP <- read.csv("../PFAS point source data/WWTP facility_details.csv")

# landfills from LMOP
landfills.lmop <- read.csv("../PFAS point source data/Secondary analysis/landfilllmopdata_clean.csv")

################################################################################
#### Aggregate sources into HUCs
################################################################################
### WWTPs ###
# add leading zero to HUC codes if they are under 8 characters
WWTP$Watershed.HUC <- ifelse(str_length(WWTP$Watershed.HUC) < 8, 
                              paste("0", WWTP$Watershed.HUC, sep=""),
                              WWTP$Watershed.HUC)

WWTP.HUCs <- WWTP %>% group_by(Watershed.HUC, Watershed.Name) %>% 
  dplyr::summarise(WWTP_count = n(), 
         WWTP_existingtotalflow_Mgal.d = sum(as.numeric(Existing.Total.Flow..Mgal.d.))) %>% 
  ungroup()

WWTP.HUCs <- WWTP.HUCs %>% dplyr::rename(
  HUC_CODE = Watershed.HUC,
  HUC_NAME = Watershed.Name
)

### Military sites ###
# note: 9 of those dropped from here are in AK
military.join <- st_intersection(military, HUCs)

# count number in each HUC
military.HUCs <- military.join %>% st_drop_geometry() %>% 
  group_by(HUC_CODE, HUC_NAME) %>%
  dplyr::summarise(MFTA_count = n()) %>% ungroup()

### EPA Stewardship sites ###
# separate into lat and lon
epa <- separate(epa, col = "Coordinates", into = c("lat", "lon"), sep = ",")

# convert to sf object
epa <- st_as_sf(epa, coords = c("lon", "lat"), crs = "NAD83")
epa <- st_transform(epa, crs = st_crs(HUCs))

# aggregate into HUCs
epa.joining <- st_intersection(epa, HUCs)

epa.HUCs <- epa.joining %>% st_drop_geometry() %>% 
  group_by(HUC_CODE, HUC_NAME) %>%
  dplyr::summarise(industries_count = n()) %>% ungroup()

### Airports ###
# clean up data
airports.trim <- airports %>% dplyr::select(LocationID, Region, State, StateName, County, 
                                CountyState, City, FacilityName, ARPLatitude, 
                                ARPLatitudeS, ARPLongitude, ARPLongitudeS)

airports.trim <- separate(airports.trim, col = "ARPLatitude", 
                     into = c("lat.deg", "lat.min", "lat.sec.dir"), sep = "-")

airports.trim <- separate(airports.trim, col = "ARPLongitude", 
                     into = c("lon.deg", "lon.min", "lon.sec.dir"), sep = "-")

# create columns without direction
airports.trim <- airports.trim %>% 
  mutate(lat.sec = as.numeric(substr(lat.sec.dir, 1, nchar(lat.sec.dir)-1)),
         lon.sec = as.numeric(substr(lon.sec.dir, 1, nchar(lon.sec.dir)-1))) 

# calculate lat and lon
airports.trim <- airports.trim %>% 
  mutate(
    lat = ifelse(grepl("N", lat.sec.dir), 
                      as.numeric(lat.deg) + as.numeric(lat.min)/60 + lat.sec/3600,
                      -1 * (as.numeric(lat.deg) + as.numeric(lat.min)/60 + lat.sec/3600)),
    lon = ifelse(grepl("E", lon.sec.dir), 
                      as.numeric(lon.deg) + as.numeric(lon.min)/60 + lon.sec/3600,
                      -1 * (as.numeric(lon.deg) + as.numeric(lon.min)/60 + lon.sec/3600)),
  )

# convert to sf object
airports.trim <- st_as_sf(airports.trim, coords = c("lon", "lat"), crs = "NAD83")
airports.trim <- st_transform(airports.trim, crs = st_crs(HUCs))

# aggregate into HUCs
# note: 47 of those dropped from here are not in the continental US
airports.joining <- st_intersection(airports.trim, HUCs)

airports.HUCs <- airports.joining %>% 
  st_drop_geometry() %>% 
  group_by(HUC_CODE, HUC_NAME) %>%
  dplyr::summarise(airport_count = n()) %>% 
  ungroup()

  # now with LMOP data
  landfills.lmop.geo <- landfills.lmop %>% 
    select(Landfill.ID, Landfill.Name, 
           wasteinplace.tons = WasteinPlace_tons, lon = Longitude,
           lat = Latitude) %>% 
    filter(!is.na(lat) & !is.na(lon)) %>% 
    mutate(wasteinplace.tons = as.numeric(gsub(",", "", wasteinplace.tons))) %>%
    # drop missing coords
    st_as_sf(coords = c("lon", "lat"), crs = "NAD83")
  
  landfills.lmop.geo <- 
    st_transform(landfills.lmop.geo, 
                 crs = st_crs(HUCs))
  
  # drops those outside of continental US
  landfills.lmop.joining <- st_intersection(landfills.lmop.geo, HUCs)
  
  landfills.lmop.HUCs <- landfills.lmop.joining %>% st_drop_geometry() %>%
    group_by(HUC_CODE, HUC_NAME) %>%
    dplyr::summarise(landfill.LMOP_count = n(),
                     totalWIP.tons = as.integer(sum(wasteinplace.tons))) %>%
    ungroup()

################################################################################
# Create spatial dataset with PFAS concentrations
################################################################################
all.measurements <- left_join(all.measurements, 
                              CWS %>% select(PWSID, HUC_NAME, HUC_CODE)
                              )

# separate out systems with multiple HUCs
all.measurements <- all.measurements %>%
  separate_rows(c(HUC_NAME, HUC_CODE), sep = ",")


HUC.PFAS <- all.measurements %>%
  group_by(HUC_NAME, HUC_CODE) %>%
  summarise(across(c(PFOA, PFOS, PFNA, PFHxS, PFBS, PFHpA, ends_with("detect")), 
                   ~max(.x, na.rm = T))
            ) %>%
  ungroup()

# create PFAS dataset of HUCs
HUC.PFAS <- left_join(HUC.PFAS, HUCs) 
HUC.PFAS <- st_as_sf(HUC.PFAS, sf_column_name = "geometry")

HUC.PFAS <- HUC.PFAS %>%
  left_join(WWTP.HUCs) %>%
  left_join(epa.HUCs) %>%
  left_join(military.HUCs) %>%
  left_join(airports.HUCs) %>%
  left_join(landfills.lmop.HUCs)

# zero out those without sources
HUC.PFAS <- HUC.PFAS %>%
  mutate(across(c(MFTA_count, airport_count, industries_count,
                  WWTP_count, WWTP_existingtotalflow_Mgal.d,
                  landfill.LMOP_count, totalWIP.tons),
                ~ifelse(is.na(.x), 0, .x))
         )

HUC.PFAS <- HUC.PFAS %>%
  mutate(across(c(MFTA_count, airport_count, industries_count),
                ~ifelse(.x > 0, "One or more", "None"),
                .names = "any_{.col}")) %>%
  rename(any_MFTA = any_MFTA_count,
         any_airport = any_airport_count,
         any_industry = any_industries_count)

# any_detect variable: this is only assessed among systems that measured all 5 at least once
HUC.PFAS <- HUC.PFAS %>% mutate(
  Any_detect = ifelse( (!is.na(PFOA_detect) & !is.na(PFOS_detect) & !is.na(PFNA_detect) & !is.na(PFHxS_detect) & !is.na(PFBS_detect)) & 
                         ((PFOA_detect > 0 | PFOS_detect > 0 | PFNA_detect > 0 | PFHxS_detect > 0 | PFBS_detect > 0)), 1, 0)
)

# NA -Inf values for PFAS (these indicate no measurements)
HUC.PFAS <- HUC.PFAS %>%
  mutate(across(c(PFOA, PFOS, PFNA, PFHxS, PFBS, PFHpA, 
                  PFNA_detect, PFBS_detect, PFHxS_detect),
                ~ifelse(is.infinite(.x), NA, .x))
         )


# log all outcome data
HUC.PFAS <- HUC.PFAS %>%
  mutate(across(c(PFOA, PFOS, PFNA, PFHxS, PFBS, PFHpA),
                log, .names = "log{.col}")
         
         )

HUC.PFAS$WWTP_logflow <- ifelse(HUC.PFAS$WWTP_existingtotalflow_Mgal.d == 0 , 0,
                               log(HUC.PFAS$WWTP_existingtotalflow_Mgal.d))

# recombine into unique HUCs only
unique.HUCs <- HUC.PFAS %>%
  group_by(HUC_CODE) %>%
  summarise(geometry = st_union(geometry),
            AREA = sum(AREA),
            PERIMETER = sum(PERIMETER)) %>%
  ungroup()

HUC.PFAS <- left_join(unique.HUCs, st_drop_geometry(HUC.PFAS) %>% select(-AREA, -PERIMETER))
HUC.PFAS <- unique(HUC.PFAS)

st_drop_geometry(HUC.PFAS) %>%
  summarise(across(c(airport_count, MFTA_count, airport_count, industries_count,
                     WWTP_count, landfill.LMOP_count),
                   ~sum(.x))
            )

# these groups refer to CWS with source water in multiple HUCs
HUC.PFAS <- subset(HUC.PFAS, !st_is_empty(geometry))

if (FALSE){
st_write(HUC.PFAS,
         delete_dsn = TRUE,
         "../Merged PFAS data/HUCs_for_analysis.geojson")
}

################################################################################
# Regression modelling
################################################################################
# first a linear reg as a test for autocorrelation (this drops nothing)
test_dat <- subset(HUC.PFAS, !is.na(PFOA) & !st_is_empty(geometry))

m1a <- lm(logPFOA ~ MFTA_count + airport_count + industries_count + WWTP_count,
         data = test_dat)

test_dat$resid <- resid(m1a)

## there is evidence of spatial autocorrelation
  nb <- poly2nb(test_dat, queen = TRUE) # queen's neighbors
  lw <- nb2listw(nb, zero.policy = TRUE)
  # evaluating spatial autocorrelation
  (morans <- moran.mc(test_dat$resid, lw, nsim = 999, zero.policy = TRUE,
                     alternative="greater")) # evidence of spatial autocorrelation; positive Moran's I, p < 0.01

  # plot(morans)
  
# so, run spatial error models; using neighbors via Queen's contiguity, weights generated above
# spatial error model
m1a <- errorsarlm(logPFOA ~ MFTA_count + airport_count + industries_count + WWTP_count,
                 listw = lw, zero.policy = T, 
                 data = test_dat)

# mark those without neighbors
HUC.PFAS$noneigh <- card(nb) == 0L

  test_dat <- subset(HUC.PFAS, noneigh == F)
  nb <- poly2nb(test_dat, queen = TRUE) # queen's neighbors
  lw <- nb2listw(nb, zero.policy = TRUE)
  m1b <- errorsarlm(logPFOA ~ MFTA_count + airport_count + industries_count + WWTP_count,
                   listw = lw, zero.policy = T, 
                   data = test_dat)

# now let's cycle through all the outcomes
PFAS.vars <- c("logPFOA", "logPFOS", "logPFNA", "logPFBS", "logPFHxS")

# this function runs and extracts model output
spatial_model.f <- function(PFAS.outcome = NULL, model.formula = NULL, exclude.NDs = NULL, exclude.zeros = NULL) {
  
  if (exclude.NDs == FALSE & exclude.zeros == FALSE) {
    
    model.data <- subset(HUC.PFAS, !is.na(eval(as.name(paste(PFAS.outcome)))))
    nb <- poly2nb(model.data)
    lw <- nb2listw(nb, zero.policy = TRUE)
    
  }
  
  else if (exclude.NDs == TRUE & exclude.zeros == FALSE) {
    
  model.data <- subset(HUC.PFAS, eval(as.name(paste(PFAS.outcome))) > log(5))
  nb <- poly2nb(model.data)
  lw <- nb2listw(nb, zero.policy = TRUE)
  
  }
  
  # else exclude observations w/o neighbors
  else {
    
    model.data <- subset(HUC.PFAS, !is.na(eval(as.name(paste(PFAS.outcome)))) &
                                   noneigh == F)
    nb <- poly2nb(model.data)
    lw <- nb2listw(nb, zero.policy = TRUE)
    
  }
  
  spatial.model <- errorsarlm(paste(PFAS.outcome, model.formula),
                              listw = lw, zero.policy = TRUE,
                              data = model.data)
  
  SE <- tibble(
    dep.var = c("lambda", names(spatial.model$rest.se)),
    SE = c(spatial.model$lambda.se, spatial.model$rest.se)
  )
  
  total.var <- as.numeric( var(st_drop_geometry(model.data[, paste(PFAS.outcome)])) )
  
  R2 <- (total.var - spatial.model$s2)/total.var
  
  # summary table
  neat.spatial.model <- tibble(
    outcome = rep(PFAS.outcome, length(coef(spatial.model))),
    outcome.name = str_replace(outcome, pattern = "log", replacement = ""),
    dep.var = names(coef(spatial.model)),
    coefficient = coef(spatial.model),
    percent.change = round((exp(coefficient) - 1)*100, 1),
    coef.se = SE$SE[mapply(grepl, dep.var, SE$dep.var)],
    UCI = round((exp(coefficient + 1.96*coef.se) - 1)*100, 1),
    LCI = round((exp(coefficient - 1.96*coef.se) - 1)*100, 1),
    t.value = coefficient/coef.se,
    # using z approximation
    p.value = ifelse(round(2*pnorm(q = abs(t.value), lower.tail = F), 3) < 0.001,
                     "<0.001", round(2*pnorm(q = abs(t.value), lower.tail = F), 3)),
    AIC = AIC(spatial.model),
    R2 = R2,
    n.obs = length(spatial.model$residuals),
    formula = paste(model.formula),
    dataset = ifelse(exclude.NDs == T, "Excludes NDs", "Includes NDs")
  )
  
}

# secondary models (with WWTP flow and landfills)
# f1 <-  "~ MFTA_count + airport_count + industries_count + WWTP_logflow"
# f3 <-  "~ MFTA_count + airport_count + industries_count + WWTP_count + landfill.LMOP_count"
f1 <-  "~ MFTA_count + airport_count + industries_count + WWTP_logflow + landfill.LMOP_count"
f2 <-  "~ MFTA_count + airport_count + industries_count + WWTP_logflow"
f3 <- "~ MFTA_count + airport_count + industries_count + WWTP_count" # identical model to Hu et al., 2016

# main model: with WWTP_logflow and landfill_count
HUC_main.models <- map_dfr(PFAS.vars, ~spatial_model.f(.x, model.formula = f1, exclude.NDs = FALSE, exclude.zeros = FALSE))
HUC_sens <- map_dfr(PFAS.vars, ~spatial_model.f(.x, model.formula = f1, exclude.NDs = FALSE, exclude.zeros = TRUE))
HUC_nolandfills <- map_dfr(PFAS.vars, ~spatial_model.f(.x, model.formula = f2, exclude.NDs = FALSE, exclude.zeros = FALSE)) # more similar to Hu et al., 2016
HUC_replicate <- map_dfr(PFAS.vars, ~spatial_model.f(.x, model.formula = f3, exclude.NDs = FALSE, exclude.zeros = FALSE)) # most similar to Hu et al., 2016

# note that excluding HUCs w/o neighbors doesn't change coefficients much
# (compare HUC_main.models to HUC_sens)

HUC_main.models <- HUC_main.models %>%
  group_by(outcome) %>%
  arrange(desc(coefficient)) %>%
  mutate(coefficient.order = ifelse(dep.var != "lambda" & dep.var != "(Intercept)",
                                    row_number(outcome), NA)) %>%
  ungroup()

HUC_nolandfills <- HUC_nolandfills %>%
  group_by(outcome) %>%
  arrange(desc(coefficient)) %>%
  mutate(coefficient.order = ifelse(dep.var != "lambda" & dep.var != "(Intercept)",
                                    row_number(outcome), NA)) %>%
  ungroup()

HUC_exclude.NDs <- map_dfr(PFAS.vars, ~spatial_model.f(.x, model.formula = f1, exclude.NDs = TRUE, exclude.zeros = FALSE))

HUC_exclude.NDs <- HUC_exclude.NDs %>%
  group_by(outcome) %>%
  arrange(desc(coefficient)) %>%
  mutate(coefficient.order = ifelse(dep.var != "lambda" & dep.var != "(Intercept)",
                                    row_number(outcome), NA)) %>%
  ungroup()

################################################################################
# model checking to evaluate if using full NDs matters for regression assumptions
# this function runs and extracts model residuals, etc
check_model.f <- function(PFAS.outcome = NULL, model.formula = NULL, exclude.NDs = NULL) {
  
  # run models again...
  if (exclude.NDs == FALSE) {
    
    model.data <- subset(HUC.PFAS, !is.na(eval(as.name(paste(PFAS.outcome)))))
    nb <- poly2nb(model.data)
    lw <- nb2listw(nb, zero.policy = TRUE)
    
  }
  
  else{  
    
    model.data <- subset(HUC.PFAS, eval(as.name(paste(PFAS.outcome))) > log(5))
    nb <- poly2nb(model.data)
    lw <- nb2listw(nb, zero.policy = TRUE)
    
  }
  
  spatial.model <- errorsarlm(paste(PFAS.outcome, model.formula),
                              listw = lw, zero.policy = TRUE,
                              data = model.data)
  
  predictions <- as.numeric(predict(spatial.model)) 
  
  # # model check
  model.check <- tibble(y = as.numeric(spatial.model$y),
                        residual = as.numeric(spatial.model$residuals),
                        outcome = paste(PFAS.outcome),
                        prediction = predictions,
                        outcome.name = paste("log(",
                          str_replace(outcome, pattern = "log", replacement = ""),
                          ")", sep = ""),
                        var = model.data$WWTP_logflow,
                        dataset = ifelse(exclude.NDs == T, "Excludes NDs", "Includes NDs"))

  return(model.check)
}

  # comparing all estimates
  all.estimates <- rbind(HUC_main.models %>% mutate(dataset = "Includes NDs (substitution)"), 
                         HUC_exclude.NDs)
  
  all.estimates$dataset <- fct_relevel(all.estimates$dataset,
                                       "Excludes NDs", "Includes NDs (substitution)")
  all.estimates <- all.estimates %>%
    mutate(dep.var.clean = case_when(dep.var == "airport_count" ~ "AFFF-certified airports",
                                     dep.var == "industries_count" ~ "Major industrial facilities",
                                     dep.var == "landfill.LMOP_count" ~ "MSW landfills",
                                     dep.var == "MFTA_count" ~ "MFTA",
                                     dep.var == "WWTP_logflow" ~ "WWTP total\nexisting effluent")
           )
  
  ggplot(all.estimates %>%
           filter(dep.var != "(Intercept)" & dep.var != "lambda" & outcome.name != "PFNA"), 
         aes(x = outcome.name, y = percent.change, group = dataset, color = dataset)) +
    geom_point(shape=21, size=2, fill="white", 
               position=position_dodge(0.6)) +
    geom_errorbar(aes(ymin=LCI, ymax=UCI), width = 0.5, size = 1, 
                  position = position_dodge(0.6)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(name = "", values = MetBrewer::met.brewer("Lakota")[4:6]) +
    facet_wrap(~dep.var.clean, scales = "free") +
    jtools::theme_nice() +
    labs(x = "", y = "% change in concentration") +
    theme(strip.text = element_text(size = 12),
          axis.text = element_text(size=10),
          axis.title = element_text(size=12),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom",
          panel.spacing = unit(2, "lines")
    )
  
  if (FALSE) {
    ggsave("../Figures/HUC check/estimates_comparison.png", dpi = 400, height = 8, width = 10, bg = "white")
    write.csv(all.estimates, "../Regressions/HUC_model_results.csv")
    write.csv(HUC_replicate, "../Regressions/HUC_model_results_forcomparison.csv")
    write.csv(HUC_replicate_ENDs, "../Regressions/HUC_model_results_forcomparison_ENDs.csv")
  }
  