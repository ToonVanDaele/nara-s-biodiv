#------------------------------------------------------------------------------
# load existing models and project with new env data
#
# Input: existing models created with biomod
# Output: projections stored in the working directory: \cr
# Parameters: modeldir = location of the modelresults (biomod)
#             modelname = name of the model run
#------------------------------------------------------------------------------

projectbiomod <- function(modelname = "test"){

  library(biomod2)
  library(tidyverse)
  library(RSQLite)
  source("11_NARA_S_Functions.R")
  source("12_NARA_S_Runbiomod.R")

  modeldir <- "../data/models"

  cat("Start - ", date(), "\n",
      "model directory: ", modeldir, "\n",
      "Modelname: ", modelname, "\n")

  #---------------------------------------------------------------------------
  # Load data
  #---------------------------------------------------------------------------
  path <- paste0("../data/data-out/", modelname, "/")

  df_proj_in <- readRDS(paste0(path, "df_proj_in.RDS"))
  explspec <- readRDS(paste0(path, "explspec.RDS"))

  db <- dbConnect(SQLite(), dbname = paste0(path, "speclist.sqlite"))
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

  origdir <- getwd()
  setwd(modeldir)

  pwalk(.l = projs, .f = function(proj, sp.n, modelname){

    cat(modelname, " - ", proj, " - ", sp.n, "\n")

    explvar <- explspec %>%
      filter(sp_n == sp.n) %>%
      dplyr::select(varname) %>%
      pull()

    projdata <- filter(df_proj_in, scen == proj)

    sink("Out_proj.txt", append = TRUE)

    project.return <- biomod_projection(sp_n = sp.n,
                                   model_name = modelname,
                                   my_proj = projdata[, explvar],
                                   my_proj_xy = projdata[, c("X", "Y")],
                                   proj_name = paste0(modelname, "_", proj))
    sink()

    })

    cat(date())
    setwd(origdir)
}
