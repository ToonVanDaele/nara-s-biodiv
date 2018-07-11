#-----------------------------------------------------------------------------
# Do biomod modelling
#-----------------------------------------------------------------------------

Biomod_modeling <- function(sp_n,
                           model_name,
                           my_resp,
                           my_expl,
                           my_resp_xy) {

  require(biomod2)
  require(tidyverse)

  # Formatting data
  my_biomod_data <- BIOMOD_FormatingData(resp.var = my_resp,
                                       expl.var = my_expl,
                                       resp.xy = my_resp_xy,
                                       resp.name = sp_n,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 500,
                                       PA.strategy = "random")

  # Set modelling options
  my_biomodoption <- BIOMOD_ModelingOptions()
  my_biomodoption <- BIOMOD_ModelingOptions(
    GAM = list(algo = "GAM_mgcv",
               myFormula = gam_formula(sp_n, names(my_expl)),
               k = -1,
               method = "GCV.Cp",
               select = TRUE,
               knots = NULL,
               paramPen = NULL))

  # Modeling
  my_biomodmodel_out <- BIOMOD_Modeling(
    my_biomod_data,
    models = c("GAM", "RF"),
    models.options = my_biomodoption,
    NbRunEval = 10,
    DataSplit = 70,
    Prevalence = 0.5,
    VarImport = 0,
    models.eval.meth = c("TSS", "ROC", "KAPPA"),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = FALSE,
    modeling.id = model_name)

  # Modeling ensemble
  my_biomod_em <- BIOMOD_EnsembleModeling(modeling.output = my_biomodmodel_out,
                           chosen.models = "all",
                           em.by = "PA_dataset+algo",
                           eval.metric = NULL,
                           eval.metric.quality.threshold = NULL,
                           models.eval.meth = c("TSS", "ROC", "KAPPA"),
                           prob.mean = TRUE,
                           prob.cv = FALSE,
                           prob.ci = FALSE,
                           prob.ci.alpha = 0.05,
                           prob.median = FALSE,
                           committee.averaging = FALSE,
                           prob.mean.weight = FALSE,
                           VarImport = 1)

  return("ok")
}

#----------------------------------------------------------------------
# Do Biomod projections present and all future scenarios
#----------------------------------------------------------------------

biomod_projection <- function(sp_n,
                             model_name,
                             my_proj,
                             my_proj_xy,
                             proj_name){
  require(tidyverse)
  require(biomod2)

  # Load the individual models
  filename <- paste0(sp_n, "/", sp_n, ".", model_name, ".models.out")
  mymodel <- loadrdata(filename)

  # Projection of individual models
  my_biomod_projection <- BIOMOD_Projection(
    modeling.output = mymodel,
    new.env = my_proj,
    proj.name = proj_name,
    xy.new.env = my_proj_xy,
    selected.models = "all",
    binary.meth = "KAPPA",
    filtered.meth = NULL,
    compress = "xz",
    build.clamping.mask = TRUE,
    silent = TRUE)

  # Projection of the ensemble
  filename <- paste0(sp_n, "/", sp_n, ".", model_name, "ensemble.models.out")
  my_em_model <- loadrdata(filename)

  my_biomod_em_forecast <- BIOMOD_EnsembleForecasting(EM.output = my_em_model,
                                projection.output = my_biomod_projection,
                                new.env = NULL,
                                xy.new.env = NULL,
                                selected.models = "all",
                                proj.name = NULL,
                                binary.meth = "KAPPA",
                                filtered.meth = NULL,
                                compress = TRUE)

}
