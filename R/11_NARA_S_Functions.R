# NARA_S_Functions
#
#------------------------------------------------------------------------
#' Get the probability and binary output from biomod
#'
#' Retrieves the probabilities and binary output
#'
#------------------------------------------------------------------------
getprobs <- function(sp_n, scenario, projname, df_proj_in){

  # load ensemble mean probability projection
  file.name <- paste0("../data/models/", sp_n, "/proj_", projname, "/proj_",
                      projname, "_", sp_n, "_ensemble.RData")

  proj.prob <- loadrdata(file.name)

  probs <- as.vector(proj.prob[, 1]) # Get the mean ensemble probability
  remove(proj.prob)

  #load ensemble binary projection
  file.name <- paste0("../data/models/", sp_n, "/proj_", projname, "/proj_",
                      projname, "_", sp_n, "_ensemble_ROCbin.RData")

  proj.bin <- loadrdata(file.name)
  bins <- as.vector(proj.bin) # Get the mean binary output of the ensemble
  remove(proj.bin)

  # bind utmID, 'probs' and 'binaries'
  idxy <- df_proj_in %>%
    filter(scen == scenario) %>%
    select(utmID)
  df.temp <- cbind(idxy, projname, sp_n, probs, bins)

  return(df.temp)
}

#------------------------------------------------------------------------
#' Get evaluation data from biomod output
#'
#' Retrieves the evaluation data from the biomod output files
#'@param sp_n species name
#'@param modelname name of the model
#'@return data frame with evaluation data
#------------------------------------------------------------------------
getevaluation <- function(sp_n, modelname) {

  require(biomod2)

    filename <- paste0("../data/models/", sp_n, "/", sp_n, ".", modelname,
                       ".models.out")
  temp <- load(filename)
  mymodel <- get(temp)
  rm(list = c(temp, "temp"))

  df <- get_evaluations(mymodel, as.data.frame = TRUE)
  df$sp_n <- sp_n
  df$modelname <- modelname

  return(df)
}

#-----------------------------------------------------------------
#' Load an RData file and return
#'
#' Loads .Rdata file
#' @param fileName file name
#-----------------------------------------------------------------
loadrdata <- function(fileName){

  load(fileName)
  get(ls()[ls() != "fileName"])
}

#------------------------------------------------------------------
# 'Create continuous maps  (png)
#
# '@param df input dataframe
# '@param path path for saving
# '@param titlecols column names from df to be used in the title
# '@param vlimits force limits c(min, max)
#------------------------------------------------------------------
plotcontinuous <- function(df, path, titlecols, vlimits) {

  require(ggplot2)
  title <- df[1, titlecols] %>%
    unite(new, titlecols, sep = "_") %>%
    pull()

  my.plot <- ggplot(df, aes(x = X, y = Y)) +
    geom_point(aes(colour = probs), size = 1) +
    scale_colour_gradientn(colours = rainbow(3), limits = vlimits) +
    coord_fixed() +
    ggtitle(title)

  ggsave(paste0(path, title, ".png"), my.plot, dpi = 150)
}

#------------------------------------------------------------------
# 'Create binary maps  (png)
#
# '@param df input dataframe
# '@param path path for saving
# '@param titlecols column names from df to be used in the title
#------------------------------------------------------------------
plotbinaries <- function(df, path, titlecols) {

  require(ggplot2)
  mycolors <- c("grey", "red")
  title <- df[1, titlecols] %>%
    unite(new, titlecols, sep = "_") %>%
    pull()

  my.plot <- ggplot(df, aes(x = X, y = Y)) +
    geom_point(aes(colour = factor(bins)), size = 1) +
    scale_color_manual(values = mycolors) +
    coord_fixed() +
    ggtitle(title)
  ggsave(paste0(path, title, ".png"), my.plot, dpi = 150)
}

#--------------------------------------------------------
# plotsum
#--------------------------------------------------------

plotsum <- function(df, path) {

  title <- df[1, "projname"]

    (p <- ggplot(df, aes(x = X, y = Y)) +
       geom_point(aes(colour = probsum), size = 1) +
       scale_colour_gradientn(colours = rainbow(3)) +
       coord_fixed()) +
       ggtitle(title)

  ggsave(paste0(path, title, "_probsum.png"), p)

  (p <- ggplot(df, aes(x = X, y = Y)) +
        geom_point(aes(colour = binsum), size = 1) +
        scale_colour_gradientn(colours = rainbow(3)) +
        coord_fixed()) +
        ggtitle(title)

  ggsave(paste0(path, title, "_binsum.png"), p)

}

#---------------------------------------------------------------------------
#' Calculate variable importance
#'
#' Calculation of the variable importance with the build in biomod function
#' @param sp_n species name
#' @param model_name model name
#' @return a data frame
#----------------------------------------------------------------------------
get_variable_importance <- function(sp_n, modelname){

  orig_wd <- getwd()
  setwd("../data/models")

  filename <- paste0(sp_n, "/", sp_n, ".", modelname, "ensemble.models.out")

  myEM <- loadrdata(filename)
  myEMvi <- get_variables_importance(myEM)

  myEMvi <- plyr::adply(myEMvi, 1:3)
  colnames(myEMvi) <- c("Varname", "rand", "EMmodel", "vi")
  myEMvi$Varname <- as.character(myEMvi$Varname)
  myEMvi$rand <- as.character(myEMvi$rand)
  myEMvi$EMmodel <- as.character(myEMvi$EMmodel)
  myEMvi$sp_n <- sp_n

  setwd(orig_wd)

  return(myEMvi)
}

#---------------------------------------------------------------------------
#' Create response plot
#'
#' Plot a 2D response plot for each variable
#' @param sp_n species name
#' @param model_name model name
#' @return a data frame
#----------------------------------------------------------------------------
create_responseplot <- function(sp_n, modelname){

  filenameEM <- paste0(sp_n, "/", sp_n, ".", modelname, "ensemble.models.out")

  myEM <- loadrdata(filenameEM)
  myEMmodels <- BIOMOD_LoadModels(myEM)

  filename <- paste0(sp_n, "/", sp_n, ".", modelname, ".models.out")
  mymod <- loadrdata(filename)
  mymodels <- BIOMOD_LoadModels(mymod)

  print(head(get_formal_data(mymod, "expl.var")))

  Plot2D <- response.plot2(models = myEMmodels,
                          Data = get_formal_data(mymod, "expl.var"),
                          show.variables = get_formal_data(mymod,
                                   "expl.var.names"),
                          do.bivariate = FALSE,
                            fixed.var.metric = "median",
                          plot = TRUE,
                          data_species = get_formal_data(mymod, "resp.var"))
}

#--------------------------------
# Select by VIF value
#--------------------------------

vif_select <- function(species, explvar, spr_plant, df_expl) {
  varnames <- explvar %>%
    filter(sp_n == species) %>%
    pull(varname)

  sphokken <- spr_plant %>%
      filter(!is.na(get(species))) %>%
      pull(utmID)

  df_vif <- df_expl %>%
    dplyr::select(-vargroup) %>%
    filter(scen == "Current" & varname %in% varnames) %>%
    spread(key = varname, value = value) %>%
    filter(utmID %in% sphokken) %>%
    dplyr::select(-utmID, -X, -Y, -scen)

  vifout <- vifstep(df_vif, th = 5)

  vifvalues <- vifout@results
  vifvalues$sp_n <- species
  vifvalues$Variables <- as.character(vifvalues$Variables)
  vifvalues <- dplyr::rename(vifvalues, varname = Variables)
  return(vifvalues)
}

#----------------------------------------------------------------------------
# Get possibly relevant variables
#----------------------------------------------------------------------------
getexplvar <- function(sp_n, spr_plant, df_expl) {

  # aantal hokken met soort sp_n
  sphokken <- spr_plant %>%
    filter(!is.na(get(sp_n))) %>%
    pull(utmID)

  # maximaal aantal variabelen (voldoende presences als absences vereist)
  maxnbvar <- trunc(min(length(sphokken),
                        nrow(spr_plant) - length(sphokken)) / 15)

  # aantal hokken met waarde > 0
  temp <- df_expl %>%
    filter(utmID %in% sphokken & !varname == "Opp_m2" & scen == "Current" &
             !value == 0) %>%
    group_by(varname) %>%
    summarise(count = n(),
              distinct = n_distinct(value))

  temp2 <- data.frame(sp_n, nbhokken = length(sphokken), maxnbvar, temp,
                      stringsAsFactors = FALSE)
}

#------------------------------------------------------------
# Create GAM formula
#------------------------------------------------------------

gam_formula <- function(spname, varnames) {

  varcs <- paste0("s(", varnames, ", k = -1, bs='cs')")

  fromstr <- paste0(spname, " ~ 1 +",
                    paste(varcs, collapse = " + "))

  form <- formula(fromstr)
  return(form)
}
