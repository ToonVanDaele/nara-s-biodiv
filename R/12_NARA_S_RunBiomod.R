#-----------------------------------------------------------------------------
# Do biomod modelling
#-----------------------------------------------------------------------------

BiomodModeling <- function(sp.n,
                           model.name,
                           my.resp,
                           my.expl,
                           my.resp.xy) {

  require(biomod2)
  require(tidyverse)

  #my.explxy <- cbind(my.expl, my.resp.xy)

  # Formatting data
  myBiomodData <- BIOMOD_FormatingData(resp.var = my.resp,
                                       expl.var = my.expl,
                                       resp.xy = my.resp.xy,
                                       resp.name = sp.n,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 500,
                                       PA.strategy = 'random')

  # Set modelling options
  myBiomodOption <- BIOMOD_ModelingOptions()
  myBiomodOption <- BIOMOD_ModelingOptions(
    GAM = list(algo = 'GAM_mgcv',
               myFormula = GAMformula(sp.n, names(my.expl)),
               k = 4,
               method = 'GCV.Cp',
               select = TRUE,
               knots = NULL,
               paramPen = NULL))

  # Modeling
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c("GAM"),
    models.options = myBiomodOption,
    NbRunEval = 10,   # Only for testing. Needs higher number
    DataSplit = 70,
    Yweights = NULL,
    VarImport = 0,
    models.eval.meth = c('TSS', 'ROC'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,      #Default is TRUE
    modeling.id = model.name)

  #myBiomodModelOut

  # Modeling ensemble
  myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                           chosen.models = 'all',
                           em.by = 'all',
                           eval.metric = 'ROC',
                           eval.metric.quality.threshold = NULL,
                           models.eval.meth = c('ROC'),
                           prob.mean = TRUE,
                           prob.cv = FALSE,
                           prob.ci = FALSE,
                           prob.ci.alpha = 0.05,
                           prob.median = FALSE,
                           committee.averaging = FALSE,
                           prob.mean.weight = FALSE,
                           prob.mean.weight.decay = 'proportional',
                           VarImport = 5)

  return("ok")
}

#----------------------------------------------------------------------
# Do Biomod projections present and all future scenarios
#----------------------------------------------------------------------

BiomodProjection <- function(sp.n,
                             model.name,
                             my.proj,
                             my.proj.xy,
                             proj.name){
  require(tidyverse)
  require(biomod2)

  # Load the individual models
  filename <- paste0(sp.n, "/", sp.n, ".", model.name, ".models.out")
  mymodel <- loadRData(filename)

  # Projection of individual models
  myBiomodProjection <- BIOMOD_Projection(
    modeling.output = mymodel,
    new.env = my.proj,
    proj.name = proj.name,
    xy.new.env = my.proj.xy,
    selected.models = 'all',
    binary.meth = 'ROC',
    filtered.meth = NULL,
    compress = 'xz',
    build.clamping.mask = TRUE,
    silent = TRUE)

  # Projection of the ensemble
  filename <- paste0(sp.n, "/", sp.n, ".", model.name, "ensemble.models.out")
  myEMmodel <- loadRData(filename)

  MyBiomodEMForecast <- BIOMOD_EnsembleForecasting(EM.output = myEMmodel,
                                projection.output = myBiomodProjection,
                                new.env = NULL,
                                xy.new.env = NULL,
                                selected.models = 'all',
                                proj.name = NULL,
                                binary.meth = 'ROC',
                                filtered.meth = NULL,
                                compress = TRUE)

}
