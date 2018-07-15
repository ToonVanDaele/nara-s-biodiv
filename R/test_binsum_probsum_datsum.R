# Compare binsum, probsum and datasum

library(tidyverse)
modelname <- "prev"
path <- paste0("../data/data-out/", modelname, "/")

df_eval <- readRDS(file = paste0(path, "df_eval.RDS"))
df_sum <- readRDS(file = paste0(path, "df_sum.RDS"))
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
  group_by(sp_n, algo, Eval.metric) %>%
  summarise(mEval = mean(Testing.data),
            mCutoff = mean(Cutoff)) %>%
  gather(key = metric, value = value, -sp_n, -Eval.metric, -algo) %>%
  unite(col = "metric", c("algo", "Eval.metric", "metric")) %>%
  spread(key = metric, value = value)

# Species prevalence
df_t2 <- df_data_in %>%
  dplyr::select(sp_n_all) %>%
  gather(key = sp_n, value = presence) %>%
  group_by(sp_n) %>%
  summarize(prevalence = sum(presence, na.rm = TRUE) )

df_tt <- left_join(df_t1, df_t2, by = "sp_n")

# ROC, TSS, KAPPA
pairs(df_tt[,c("GAM_TSS_mEval", "GAM_ROC_mEval", "GAM_KAPPA_mEval")])

# AUC/ROC and TSS versus prevalence
ggplot(df_tt, aes(x = prevalence, y = GAM_ROC_mEval)) + geom_point()
ggplot(df_tt, aes(x = prevalence, y = GAM_TSS_mEval)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = GAM_KAPPA_mEval)) + geom_point() + geom_smooth()

# Threshold (cutoff) versus presences
ggplot(df_tt, aes(x = prevalence, y = GAM_ROC_mCutoff)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = GAM_TSS_mCutoff)) + geom_point() + geom_smooth()
ggplot(df_tt, aes(x = prevalence, y = GAM_KAPPA_mCutoff)) + geom_point() + geom_smooth()

# Verify Binsum (sum of binary results), probsum (sum of probabilities),
# datasum (sum of observation data)

df_plant_data <- df_data_in %>%
  dplyr::select(utmID, X, Y, sp_n_all) %>%
  gather(key = sp_n, value = presence, -utmID, -X, -Y)
df_plant_data[is.na(df_plant_data)] <- 0

# Join evaluation and prevalence information
df_plant_data <- left_join(x = df_plant_data,
                           y = df_tt %>%
                             dplyr::select(sp_n, prevalence,
                                           GAM_ROC_mEval, GAM_KAPPA_mEval,
                                           GAM_ROC_mCutoff),
                           by = "sp_n")

df_probbin <- df_probs %>%
  filter(projname == "rlw_Current" &
           utmID %in% utmIDs)

df_all <- inner_join(x = df_plant_data,
                     y = df_probbin,
                     by = c("utmID", "sp_n"))

head(df_all)

df_all_sums <- df_all %>%
  filter(GAM_ROC_mEval >= 0.7) %>%
  group_by(utmID) %>%
  summarise(datasum = sum(presence),
            probsum = sum(probs) / 1000,
            binsum = sum(bins))

head(df_all_sums)
summary(df_all_sums)

df_all_sums %>%
  gather(key = data, value = nb, -utmID, -datasum) %>%
  ggplot(aes(x = datasum, y = nb)) +
  geom_jitter(size = 0.2) + geom_smooth(method = 'lm') +
  geom_abline(slope = 1) +
  coord_cartesian(xlim = c(0,250), ylim = c(0, 250)) +
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

df_sp_prev$lroc <- cut(x = df_sp_prev$GAM_ROC_mEval, breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))

df_sp_prev$lthres <- cut(x = df_sp_prev$GAM_ROC_mCutoff,
                         breaks = quantile(x = df_sp_prev$GAM_ROC_mCutoff))

ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs, colour = lthres)) + geom_point()

df_sp_prev %>%
  filter(GAM_ROC_mEval >= 0.7) %>%
  ggplot(aes(x = prev_data, y = prev_probs)) + geom_point() +
    geom_abline(slope = 1) + geom_smooth()

df_sp_prev %>%
  filter(GAM_ROC_mEval >= 0.7) %>%
  ggplot(aes(x = prev_data, y = prev_bins)) + geom_point() +
    geom_abline(slope = 1) + geom_smooth()

# Zoom in to 4 species on the break prevalence = 1250
ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs)) + geom_point() +
  coord_cartesian(xlim = c(1150, 1350), ylim = c(1250, 1500))

ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs, label = sp_n)) + geom_point() +
  coord_cartesian(xlim = c(1200, 1300), ylim = c(1250, 1600)) + geom_text()

df_sp_prev %>%
  filter(prev_data > 1240 & prev_data < 1295)


#

ggplot(df_sp_prev, aes(x = prev_data, y = GAM_ROC_mCutoff, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_data, y = GAM_KAPPA_mCutoff, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = prev_probs, y = prev_bins, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = prev_probs, y = GAM_ROC_mCutoff, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_bins, y = GAM_ROC_mCutoff, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = prev_bins, y = GAM_KAPPA_mCutoff, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = GAM_KAPPA_mEval, y = GAM_KAPPA_mCutoff)) + geom_point()

ggplot(df_sp_prev, aes(x = GAM_ROC_mEval, y = prev_probs - prev_data)) + geom_point()
ggplot(df_sp_prev, aes(x = GAM_ROC_mEval, y = prev_bins - prev_data)) + geom_point()

ggplot(df_sp_prev, aes(x = GAM_ROC_mEval, y = prev_data)) + geom_point()




head(df_plant_data)

ggplot()
