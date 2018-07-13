#------------------------------------------------------------------------------
# load preprocessed data, run biomod models and write output
#
# Input: df_expl_in.RDS, df_plant_in.RDS, explspec.RDS \cr
# Output: Biomod stores model output in the working directory: \cr
# Parameters: modeldir = location of the data (input, output and modelresults)
#             modelname = name of the model run
#------------------------------------------------------------------------------

procesbiomod <- function(modelname = "test") {

  library(biomod2)
  library(tidyverse)
  library(RSQLite)
  source("11_NARA_S_Functions.R")
  source("12_NARA_S_Runbiomod.R")

  modeldir <- "../data/models"
  path <- paste0("../data/data-out/", modelname, "/")

  cat("Start - ", date(), "\n",
      "Model directory: ", modeldir, "\n",
      "Modelname: ", modelname, "\n",
      "data from dir: ", path, "\n")

  #---------------------------------------------------------------------------
  # Load data
  #---------------------------------------------------------------------------

  df_data_in <- readRDS(paste0(path, "df_data_in.RDS"))
  explspec <- readRDS(paste0(path, "explspec.RDS"))

  db <- dbConnect(SQLite(), dbname = paste0(path, "speclist.sqlite"))
  speclist <- dbReadTable(db, "speclist")

  #---------------------------------------------------------------------------
  # Run biomod in a loop for all species in speclist
  # models are stored on hard disk in the working directory
  #---------------------------------------------------------------------------

  #Select only species that didn't run yet (modelrun == "-")
  my.spec <- speclist %>%
    filter(modelrun == "-") %>%
    dplyr::select("species") %>%
    pull()

  if (length(my.spec) == 0) {
    print("No models to run - check modelname")
    }

  origdir <- getwd()
  setwd(modeldir)

  walk2(.x = my.spec, .y = modelname, .f = function(sp.n, modelname){

    cat(modelname, " - ", sp.n, "\n")

    my.resp <- df_data_in[, sp.n]
    explvar <- explspec %>%
      filter(sp_n == sp.n) %>%
      dplyr::select(varname) %>%
      pull() %>%
      intersect(colnames(df_data_in))

   sink("Out.txt", append = TRUE)
    model.return <- Biomod_modeling(sp_n = sp.n,
                                   model_name = modelname,
                                   my_resp = my.resp,
                                   my_expl = df_data_in[, explvar],
                                   my_resp_xy = df_data_in[, c("X", "Y")])

    sink()

    dbSendQuery(db, paste0("UPDATE speclist SET modelrun = '", model.return,
                     "' WHERE species = '", sp.n, "'"))
    })

    dbDisconnect(db)
    cat(date())
    setwd(origdir)
}
