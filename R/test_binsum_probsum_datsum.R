# Compare binsum, probsum and datasum


library(tidyverse)
modelname <- "allsp"
path <- paste0("../data/data-out/", modelname, "/")

df_eval <- readRDS(file = paste0(path, "df_eval.RDS"))
df_sums_AUC <- readRDS(file = paste0(path, "df_sums_AUC.RDS"))
df_data_in <- readRDS(file = paste0(path, "df_data_in.RDS"))
df_probs <- readRDS(file = paste0(path, "df_probs.RDS"))

head(df_probs)

temp <- df_probs %>%
  group_by(utmID, projname) %>%
  summarise(bin = sum(bins),
            prob = sum(probs) / 1000)

summary(temp)

temp2 <-  filter(temp, projname == "allsp_Current") %>%
  dplyr::select(utmID, curbin = bin, curprob = prob)

temp2$probcut <- cut(temp2$curprob, breaks = c(0, 243, 272, 298, 385))
temp2$bincut <- cut(temp2$curbin, breaks = c(0, 193, 245, 297, 466))

summary(temp2)

temp <- left_join(x = temp,
                  y = temp2 %>%
                    select(probcut, bincut, utmID, curbin, curprob),
                  by = "utmID")

summary(temp)

temp$binverschil <- temp$bin - temp$curbin
temp$probverschil <- temp$prob - temp$curprob

temp %>%
  filter(!projname == "allsp_Current") %>%
  ggplot(aes(x = projname, y = probverschil)) + geom_boxplot() +
  geom_hline(yintercept = 0) +
    facet_wrap(~probcut)

temp %>%
  ggplot(aes(x = projname, y = prob)) + geom_boxplot() +
    facet_wrap(~probcut)

temp %>%
  ggplot(aes(x = projname, y = prob)) + geom_boxplot()

temp %>%
  filter(!projname == "allsp_Current") %>%
  ggplot(aes(x = projname, y = probverschil)) + geom_boxplot() +
  geom_hline(yintercept = 0)

ggplot(temp, aes(x = prob, colour = projname)) + geom_density()


# total list of species
sp_n_all <- unique(df_eval$sp_n)

# list of utmID with data
utmIDs <- unique(df_data_in$utmID)

# AUC versus presences for all models

# Evaluation metrics
df_t1 <- df_eval %>%
  filter(Eval.metric == "ROC") %>%
  group_by(sp_n) %>%
  summarise(mROC = mean(Testing.data),
            mCutoff = mean(Cutoff),
            mSensitivity = mean(Sensitivity),
            mSpecificity = mean(Specificity))

# Species prevalence
df_t2 <- df_data_in %>%
  select(sp_n_all) %>%
  gather(key = sp_n, value = presence) %>%
  group_by(sp_n) %>%
  summarize(count = sum(presence, na.rm = TRUE))

df_tt <- left_join(df_t1, df_t2, by = "sp_n")

# AUC versus prevalence
ggplot(df_tt, aes(x = count, y = mROC)) + geom_point() + geom_smooth()

df_tt %>%
  filter(mROC >=0.7) %>%
  ggplot(aes(x = count, y = mROC)) + geom_point() + geom_smooth()

# Threshold (cutoff) versus presences
ggplot(df_tt, aes(x = count, y = mCutoff)) + geom_point() + geom_smooth()

# Senisitity versus presences
ggplot(df_tt, aes(x = count, y = mSensitivity)) + geom_point() + geom_smooth()

# Specificity versus presences
ggplot(df_tt, aes(x = count, y = mSpecificity)) + geom_point() + geom_smooth()


# Verify Binsum (sum of binary results), probsum (sum of probabilitie),
# datasum (sum of observation data)

df_plant_data <- df_data_in %>%
  select(utmID, X, Y, sp_n_all) %>%
  gather(key = sp_n, value = presence, -utmID, -X, -Y)
df_plant_data[is.na(df_plant_data)] <- 0

df_probbin <- df_probs %>%
  filter(projname == "allsp_Current" &
           utmID %in% utmIDs)

df_all <- inner_join(x = df_plant_data,
                     y = df_probbin,
                     by = c("utmID", "sp_n"))

head(df_all)

df_all_sums <- df_all %>%
  filter(mROC >= 0.7) %>%
  group_by(utmID) %>%
  summarise(datasum = sum(presence),
            probsum = sum(probs) / 1000,
            binsum = sum(bins))

head(df_all_sums)

df_all_sums %>%
  gather(key = data, value = nb, -utmID, -datasum) %>%
  ggplot(aes(x = datasum, y = nb)) + geom_point(size = 0.2) + geom_smooth(method = 'lm') +
  geom_abline(slope = 1) +
  coord_cartesian(xlim = c(0,270), ylim = c(0, 270)) +
  facet_wrap(~data)

cor(x = df_all_sums[,2:4])
cor(df_all_sums$datasum, df_all_sums$binsum)
cor(df_all_sums$datasum, df_all_sums$probsum)
cor(df_all_sums$binsum, df_all_sums$probsum)

lm1 <- lm(formula = binsum ~ datasum, data = df_all_sums)
summary(lm1)
anova(lm1)

# species count (prevalence) versus

df_sp_prev <- df_all %>%
  group_by(sp_n) %>%
  summarise(prev_data = sum(presence),
            prev_probs = sum(probs / 1000),
            prev_bins = sum(bins)) %>%
  left_join(y = df_tt,
            by = "sp_n")

df_sp_prev$lroc <- cut(x = df_sp_prev$mROC, breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1))


ggplot(df_sp_prev, aes(x = prev_data, y = prev_probs, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = prev_data, y = prev_bins, colour = lroc)) + geom_point()

ggplot(df_sp_prev, aes(x = mROC, y = prev_probs - prev_data, colour = lroc)) + geom_point()
ggplot(df_sp_prev, aes(x = mROC, y = prev_bins - prev_data)) + geom_point()

