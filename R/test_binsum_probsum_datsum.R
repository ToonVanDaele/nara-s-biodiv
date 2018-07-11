# Compare binsum, probsum and datasum


library(tidyverse)
modelname <- "kappa"
path <- paste0("../data/data-out/", modelname, "/")

df_eval <- readRDS(file = paste0(path, "df_eval.RDS"))
df_sums_AUC <- readRDS(file = paste0(path, "df_sums_AUC.RDS"))
df_data_in <- readRDS(file = paste0(path, "df_data_in.RDS"))
df_probs <- readRDS(file = paste0(path, "df_probs.RDS"))

head(df_probs)

# total list of modelled species
sp_n_all <- unique(df_eval$sp_n)

# list of utmID with plant data
utmIDs <- unique(df_data_in$utmID)
nbutm <- length(utmIDs)

# Evaluation metrics per species model TSS, ROC, Kappa
df_t1 <- df_eval %>%
  group_by(sp_n, Eval.metric) %>%
  summarise(mEval = mean(Testing.data),
            mCutoff = mean(Cutoff),
            mSens = mean(Sensitivity),
            mSpec = mean(Specificity)) %>%
  gather(key = metric, value = value, -sp_n, -Eval.metric) %>%
  unite(col = "metric", c("metric", "Eval.metric")) %>%
  spread(key = metric, value = value)

# Species prevalence
df_t2 <- df_data_in %>%
  dplyr::select(sp_n_all) %>%
  gather(key = sp_n, value = presence) %>%
  group_by(sp_n) %>%
  summarize(prevalence = sum(presence, na.rm = TRUE) / nbutm)

df_tt <- left_join(df_t1, df_t2, by = "sp_n")

# ROC, TSS, KAPPA
pairs(df_tt[,c("mEval_TSS", "mEval_ROC", "mEval_KAPPA")])

# AUC/ROC and TSS versus prevalence
ggplot(df_tt, aes(x = prevalence, y = mEval_ROC)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = mEval_TSS)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = mEval_KAPPA)) + geom_point() + geom_smooth()

# Threshold (cutoff) versus presences
ggplot(df_tt, aes(x = prevalence, y = mCutoff_ROC)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = mCutoff_TSS)) + geom_point() + geom_smooth()

# Senisitity versus presences
ggplot(df_tt, aes(x = prevalence, y = mSens_ROC)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = mSens_TSS)) + geom_point() + geom_smooth()

# Specificity versus presences
ggplot(df_tt, aes(x = prevalence, y = mSpec_ROC)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = mSpec_TSS)) + geom_point() + geom_smooth()

# Verify Binsum (sum of binary results), probsum (sum of probabilitie),
# datasum (sum of observation data)

df_plant_data <- df_data_in %>%
  dplyr::select(utmID, X, Y, sp_n_all) %>%
  gather(key = sp_n, value = presence, -utmID, -X, -Y)
df_plant_data[is.na(df_plant_data)] <- 0

# Join mROC and prevalence information
df_plant_data <- left_join(x = df_plant_data,
                           y = df_tt %>%
                             dplyr::select(sp_n, mEval_ROC, mEval_TSS,
                                           mEval_KAPPA, prevalence),
                           by = "sp_n")

df_probbin <- df_probs %>%
  filter(projname == "kappa_Current" &
           utmID %in% utmIDs)

df_all <- inner_join(x = df_plant_data,
                     y = df_probbin,
                     by = c("utmID", "sp_n"))

head(df_all)

df_all_sums <- df_all %>%
  filter(mEval_ROC >= 0.7) %>%
  group_by(utmID) %>%
  summarise(datasum = sum(presence),
            probsum = sum(probs) / 1000,
            binsum = sum(bins))

head(df_all_sums)

df_all_sums %>%
  gather(key = data, value = nb, -utmID, -datasum) %>%
  ggplot(aes(x = datasum, y = nb)) +
  geom_jitter(size = 0.2) + geom_smooth(method = 'lm') +
  geom_abline(slope = 1) +
  coord_cartesian(xlim = c(0,20), ylim = c(0, 20)) +
  facet_wrap(~data)

cor(x = df_all_sums[,2:4])

# species count (prevalence) versus

df_sp_prev <- df_all %>%
  group_by(sp_n) %>%
  summarise(prev_data = sum(presence),
            prev_probs = sum(probs / 1000),
            prev_bins = sum(bins)) %>%
  left_join(y = df_tt,
            by = "sp_n")

df_sp_prev$lroc <- cut(x = df_sp_prev$mEval_ROC, breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))

df_sp_prev$lthres <- cut(x = df_sp_prev$mCutoff_ROC,
                         breaks = quantile(x = df_sp_prev$mCutoff_ROC))

ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs, colour = lthres)) + geom_point()

df_sp_prev %>%
  #filter(mEval_ROC >= 0.7) %>%
  ggplot(aes(x = prev_data, y = prev_probs)) + geom_point() +
    geom_abline(slope = 1) + geom_smooth()

df_sp_prev %>%
  #filter(mEval_ROC >= 0.7) %>%
  ggplot(aes(x = prev_data, y = prev_bins)) + geom_point() +
    geom_abline(slope = 1) + geom_smooth()

ggplot(df_sp_prev, aes(x = prev_data, y = mCutoff_ROC, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_data, y = mCutoff_TSS, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_data, y = mCutoff_KAPPA, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = prev_probs, y = prev_bins, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = prev_probs, y = mCutoff_ROC, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_bins, y = mCutoff_ROC, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_bins, y = mCutoff_KAPPA, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = mEval_KAPPA, y = mCutoff_KAPPA)) + geom_point()

ggplot(df_sp_prev, aes(x = mEval_ROC, y = prev_probs - prev_data)) + geom_point()
ggplot(df_sp_prev, aes(x = mEval_ROC, y = prev_bins - prev_data)) + geom_point()



ggplot(df_sp_prev, aes(x = mEval_ROC, y = prev_data)) + geom_point()
ggplot(df_sp_prev, aes(x = mEval_TSS, y = prev_data)) + geom_point()


