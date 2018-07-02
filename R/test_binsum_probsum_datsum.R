# Compare binsum, probsum and datasum
#

library(tidyverse)
modelname <- "allsp"
path <- paste0("../data/data-out/", modelname, "/")

df_eval <- readRDS(file = paste0(path, "df_eval.RDS"))
df_sums_AUC <- readRDS(file = paste0(path, "df_sums_AUC.RDS"))
df_data_in <- readRDS(file = paste0(path, "df_data_in.RDS"))
df_probs <- readRDS(file = paste0(path, "df_probs.RDS"))

# Check treshold versus presences



# AUC versus presences

df_t1 <- df_eval %>%
  filter(Eval.metric == "ROC") %>%
  group_by(sp_n) %>%
  summarise(mROC = mean(Testing.data),
            mCutoff = mean(Cutoff))

df_t2 <- df_data_in[,23:658] %>%
  gather(key = sp_n, value = presence) %>%
  group_by(sp_n) %>%
  summarize(count = sum(presence, na.rm = TRUE))

df_tt <- left_join(df_t1, df_t2, by = "sp_n")

ggplot(df_tt, aes(x = count, y = mROC)) + geom_point() + geom_smooth()

sp_model_ok <- filter(df_tt, mROC >= 0.7) %>%
  .$sp_n


# Check Binsum, probsum, input data sum

df_plants <- df_data_in %>%
  select(utmID, sp_model_ok) %>%
  gather(key = sp_n, value = presence, -utmID) %>%
  filter(presence == 1 & sp_n %in% sp_model_ok) %>%
  group_by(utmID) %>%
  summarise(datasum = sum(presence))

df_probbin <- df_probs %>%
  filter(projname == "allsp_Current" &
           sp_n %in% sp_model_ok &
           utmID %in% df_plants$utmID) %>%
  group_by(utmID) %>%
  summarise(probsum = sum(probs) / 1000,
            binsum = sum(bins))

df_all <- inner_join(x = df_plants,
                     y = df_probbin,
                     by = "utmID")

df_all <- df_all %>%
  gather(key = )

ggplot(df_all, aes(x = probsum, y = binsum)) + geom_point() + geom_smooth() +
  coord_cartesian(xlim = c(20, 270), ylim = c(20, 270))

ggplot(df_all, aes(x = datasum, y = binsum)) + geom_point() + geom_smooth() +
  geom_point(aes(y = probsum), colour = "red")



