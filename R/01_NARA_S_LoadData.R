#----------------------------------------------------------------------------
# Step 1: Import plant & environmental data
#
# Loads the data from a location (absolute path);
# Some rearrangements and changes
# save data as RDS variables in '../data/data-in'
#----------------------------------------------------------------------------

options(stringsAsFactors = FALSE)
library(RODBC)
library(tidyverse)

path <- paste0("C:/Users/toon_vandaele/toon.vandaele@inbo.be/Project_WET/",
          "NARA/NARA_S/DataIN/")

#### Plant data ####

## Get data
db_plant <- odbcConnectAccess2007(access.file = paste0(path,
                                          "Nara_planten_dataset.accdb"))

df_plant <- sqlQuery(db_plant, "SELECT * FROM SoortenNaraNa2000",
                     stringsAsFactors = FALSE)
odbcClose(db_plant)

## Some checks
str(df_plant)
summary(df_plant)
head(df_plant)

## Change column names
colnames(df_plant)[1] <- "utmID"

# Replace underscore & hyphen by a dot
#(biomod doesn't like underscores or hyphen in species name)
df_plant$Taxoncode <- gsub("_", ".", df_plant$Taxoncode)
df_plant$Taxoncode <- gsub("-", ".", df_plant$Taxoncode)

#### Environmental data ####

# Current
df_cur <- read.csv2(paste0(path, "expl_var.txt"))

df_cur$scen <- "Current"
df_cur <- select(df_cur, -Shape_Length, -Shape_Area)

# KR1
df_kr1 <- read.csv2(paste0(path, "hotspotplanten_variabelen_KR1_2050.csv"))
df_kr1$scen <- "kr1"

# KR2
df_kr2 <- read.csv2(paste0(path, "hotspotplanten_variabelen_KR2_2050.csv"))
df_kr2$scen <- "kr2"

#KR3
df_kr3 <- read.csv2(paste0(path, "hotspotplanten_variabelen_KR3_2050.csv"))
df_kr3$scen <- "kr3"

#KR4
df_kr4 <- read.csv2(paste0(path, "hotspotplanten_variabelen_KR4_2050.csv"))
df_kr4$scen <- "kr4"

# Define groups of explanatory variables:
# land use (lu), soil (soil), condition (cond)
luvar <- c("bebouwing", "wegen", "recreatie", "stadsgroen", "akker",
        "prodgrasland", "heide", "loofbos", "naaldbos", "halfnatgras", "moeras",
        "duinen", "slikschor", "ovlaaggroen", "water")
soilvar <- c("textzand", "textleem", "textklei", "textveen", "textmergkrijt")
condvar <- c("drainage", "potnatzilt", "potnattrofie", "potnatzuur")

df_vargroup <- bind_rows(data.frame(vargroup = "lu", varname = luvar),
                       data.frame(vargroup = "soil", varname = soilvar),
                       data.frame(vargroup = "cond", varname = condvar))

# Some columns are missing in kr1..4 compared to the current scenario
df_vars <- data.frame(varnames = names(df_cur),
                      kr1 = names(df_cur) %in% names(df_kr1),
                      kr2 = names(df_cur) %in% names(df_kr2),
                      kr3 = names(df_cur) %in% names(df_kr3),
                      kr4 = names(df_cur) %in% names(df_kr4))

df_vars

# 'Stadsgroen' is missing in kr2 and kr3 (ok, not used in analysis)
# 'Naaldbos' is missing in kr2 and kr4 (correct, will become 0)
# 'wegen' is absent in all kr (ok, not used in analysis)
df_kr2$stadsgroen <- df_kr3$stadsgroen <- NA
df_kr2$naaldbos <- df_kr4$naaldbos <- NA
df_kr1$wegen <- df_kr2$wegen <- df_kr3$wegen <- df_kr4$wegen <- NA

# Add soil and cond data from current to kr1..kr4
varsel <- c("TAG", "POINT_X", "POINT_Y", "Opp_m2", soilvar, condvar)
df_kr1 <- left_join(df_kr1, df_cur[, varsel], by = "TAG")
df_kr2 <- left_join(df_kr2, df_cur[, varsel], by = "TAG")
df_kr3 <- left_join(df_kr3, df_cur[, varsel], by = "TAG")
df_kr4 <- left_join(df_kr4, df_cur[, varsel], by = "TAG")

df_expl <- bind_rows(df_cur, df_kr1, df_kr2, df_kr3, df_kr4)

# Drop objectID
df_expl <- select(df_expl, -OBJECTID)

# Change some column names
colnames(df_expl)[c(1, 3, 4)] <- c("utmID", "X", "Y")

# NA -> 0 for almost all variables, except 'zuurgraad' & 'trofie'
df_expl <- df_expl %>%
  mutate_at(.vars = vars(-potnatzuur, -potnattrofie),
            .funs = funs(replace(., is.na(.), 0)))

#Using m2 units give hughe numbers and numerical problems for modelling.
#Transform m2 to km2 gives values between 0 and 1
ckm <- c("bebouwing", "wegen", "recreatie", "stadsgroen", "akker",
        "prodgrasland", "heide", "loofbos", "naaldbos", "halfnatgras", "moeras",
        "duinen", "slikschor", "ovlaaggroen", "water", "potnatzilt")

df_expl[, ckm] <- df_expl[, ckm] / 1000000

# Gather data to long format
df_expl <- df_expl %>%
  gather(key = "varname", value = "value", -utmID, -X, -Y, -scen)

# Add vargroup
df_expl <- left_join(df_expl,
                      df_vargroup,
                      by = "varname")

##### Load red list species
df_redlist <- read.csv(file = paste0(path, "RL_planten.csv"), sep = ";")

#### load ecoregions utmID
df_ecoreg <- foreign::read.dbf(file = paste0(path, "ifblecoregio.dbf"),
                              as.is = TRUE)
df_ecoreg <- select(df_ecoreg, regio = REGIO, utmID = TAG)

# Save results as binary files
saveRDS(df_plant, "../data/data-in/df_plant.RDS")
saveRDS(df_expl, "../data/data-in/df_expl.RDS")
saveRDS(df_redlist, "../data/data-in/df_redlist.RDS")
saveRDS(df_ecoreg, "../data/data-in/df_ecoreg.RDS")
