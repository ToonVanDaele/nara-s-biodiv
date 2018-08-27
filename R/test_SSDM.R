#### TEST SSDM package

library(SSDM)
library(tidyverse)

path = paste0("../data/data-out/rlw/")

## ----Environmental variables---------------------------------------------
Env_table <- readRDS(paste0(path, "df_proj_in.RDS"))
Env_table <- Env_table %>% filter(scen == "Current")

myras <- raster::raster(xmn = min(Env_table$X) - 500,
               xmx = max(Env_table$X) + 500,
               ymn = min(Env_table$Y) - 500,
               ymx = max(Env_table$Y) + 500,
               resolution = 1000)

lnames <- colnames(Env_table)[5:22]

rs <- raster::stack()
for (lname in lnames){
  rl <- raster::rasterize(Env_table[,2:3], myras, field = Env_table[,lname], fun = mean)
  rs <- raster::stack(rs, rl)
}

names(rs) <- lnames
rs


## ----Occurences---------------------------------------------

df_plant <- readRDS("../data/data-in/df_plant_eda.RDS")

df_plant2 <- df_plant %>%
  left_join(dplyr::select(Env_table,
                          utmID, X, Y),
            by = "utmID")


## Set species and selection of variables

explspec <- readRDS(paste0(path, "explspec.RDS"))

df_plant3 <- df_plant2 %>%
  filter(Taxoncode == "acer.cam")


speclist <- c("acer.cam", "calluvul", "alnusglu", "equisarv")

df_plant4 <- df_plant2 %>%
  filter(Taxoncode %in% speclist)

explvar <- explspec %>%
      filter(sp_n == spec) %>%
      dplyr::select(varname) %>%
      pull()

rs_env <- raster::subset(rs, explvar)


## ----SDM-----------------------------------------------------------------

SDM <- modelling(algorithm = 'GLM',
                 Occurrences = df_plant3, Xcol = 'X', Ycol = 'Y',
                 Env = rs_env,
                 verbose = TRUE)
plot(SDM@projection, main = 'SDM - acer.cam - GAM algorithm')

plot(SDM@binary, main = 'SDM - acer.cam - GAM algorithm')

SDM@evaluation
SDM@variable.importance

## ----ESDM----------------------------------------------------------------
ESDM <- ensemble_modelling(algorithms = c('GLM', 'MARS'),
                           Occurrences = df_plant3, Xcol = 'X', Ycol = 'Y',
                           Env = rs_env, rep = 1,
                           ensemble.thresh = 0.5,
                           verbose = TRUE)

plot(ESDM@projection, main = 'ESDM')

ESDM@evaluation
ESDM@variable.importance

## ----SSDM----------------------------------------------------------------
SSDM <- stack_modelling(algorithms = c('GLM'),
                        Occurrences = df_plant4, Xcol = 'X', Ycol = 'Y',
                        Env = rs_env,
                        rep = 1,
                        Spcol = 'Taxoncode',
                        method = "bSSDM",
                        ensemble.thresh = c(0.5),
                        verbose = TRUE)

plot(SSDM@diversity.map, main = 'SSDM')

plot(SSDM@uncertainty)

SSDM@parameters

SSDM@algorithm.evaluation

SSDM@enms$acer.cam.Ensemble.SDM@evaluation

SSDM@evaluation

## ----SDM evaluation------------------------------------------------------
knitr::kable(ESDM@evaluation)

## ----SSDM evaluation-----------------------------------------------------
knitr::kable(SSDM@evaluation)

## ----SSDM variable importance--------------------------------------------
knitr::kable(SSDM@variable.importance)

## ----SSDM endemism-------------------------------------------------------
plot(SSDM@endemism.map, main = 'SSDM')

## ----plot----------------------------------------------------------------
# plot(SSDM)
