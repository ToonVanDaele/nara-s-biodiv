#------------------------------------------------------------------------------
# load existing models and project with new env data
#
# Input: existing models created with biomod
# Output: projections stored in the working directory: \cr
# Parameters: modeldir = location of the modelresults (biomod)
#             modelname = name of the model run
#------------------------------------------------------------------------------

library(biomod2)
library(tidyverse)
library(RSQLite)
source("11_NARA_S_Functions.R")
source("12_NARA_S_Runbiomod.R")

modeldir <- "../data/models"
modelname <- "allsp"

cat("Start - ", date(), "\n",
    "model directory: ", modeldir, "\n",
    "Modelname: ", modelname)

#---------------------------------------------------------------------------
# Load data
#---------------------------------------------------------------------------
pathin <- paste0("../data/data-in/", modelname, "/")

df_proj_in <- readRDS(paste0(pathin, "df_proj_in.RDS"))
explspec <- readRDS(paste0(pathin, "explspec.RDS"))

db <- dbConnect(SQLite(), dbname = paste0(pathin, "speclist.sqlite"))
speclist <- dbReadTable(db, "speclist")
dbDisconnect(db)

# List of projections to be made
projs <- expand.grid(proj = unique(df_proj_in$scen),
                     sp.n = speclist$species,
                     modelname = modelname,
                     stringsAsFactors = FALSE)

#---------------------------------------------------------------------------
# Run biomod projection in a loop for all species and scenarios
# The projections are stored on hard disk in the working directory
#---------------------------------------------------------------------------

setwd(modeldir)

pwalk(.l = projs, .f = function(proj, sp.n, modelname){

  cat(modelname, " - ", proj, " - ", sp.n, "\n")

  explvar <- explspec %>%
    filter(sp_n == sp.n) %>%
    dplyr::select(varname) %>%
    pull()

  projdata <- filter(df_proj_in, scen == proj)

  sink("Out_proj.txt", append = TRUE)

  project.return <- biomod_projection(sp.n = sp.n,
                                 model.name = modelname,
                                 my.proj = projdata[, explvar],
                                 my.proj.xy = projdata[, c("X", "Y")],
                                 proj.name = paste0(modelname, "_", proj))
  sink()

  })

  cat(date())
