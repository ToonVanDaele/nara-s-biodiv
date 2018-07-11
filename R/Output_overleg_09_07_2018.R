# Output overleg 09/07/2018

# Hieronder het script voor de output die we daarnet kregen

# Het path voor het inlezen van de data moet aangepast worden.

# Het best kan je de bestanden lokaal in een tijdelijke map
# op je laptop zetten.


library(tidyverse)
modelname <- "allsp"
path <- paste0("../data/data-out/", modelname, "/")

# path aan te passen naar bv. "c:/temp/"

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


