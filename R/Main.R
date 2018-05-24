# Main script running R and markdown scripts

modelname <- "test3"
pathout = paste0("../data/data-out/", modelname)

# 1 Load data
source("01_NARA_S_LoadData.R")

# 2 EDA
rmarkdown::render("02_NARA_S_EDA.Rmd", output_dir = "../data/data-in/")

# 3 Preprocessing
rmarkdown::render("03_NARA_S_PreProcessing.Rmd", output_dir = pathout,
                  params = list(modelname = modelname))

# 4 Biomod modelling
source("04_NARA_S_ProcesBiomod.R")

# 5 Post process biomod modelling
rmarkdown::render("05_NARA_S_PostProcesModels.Rmd", output_dir = pathout,
                  params = list(modelname = modelname,
                                calcvarimp = FALSE))

# 6 Biomod projections
source("06_NARA_S_ProjectBiomod.R")

# 7 Post proces projections
rmarkdown::render("07_NARA_S_PostProcesProject.Rmd", output_dir = pathout,
                  params = list(modelname = modelname))

# 8 Generate results
rmarkdown::render("08_NARA_S_Output.Rmd", output_dir = pathout,
                  params = list(modelname = modelname))

