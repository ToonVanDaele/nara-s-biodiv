## probmaps according prevalence

library(tidyverse)
modelname <- "allsp"
path <- paste0("../data/data-out/", modelname, "/")
source("11_NARA_S_Functions.R")

# load data
df_data_in <- readRDS(file = paste0(path, "df_data_in.RDS"))
df_probs <- readRDS(file = paste0(path, "df_probs.RDS"))
df_proj_in <- readRDS(paste0(path, "df_proj_in.RDS"))
df_redlist <- readRDS("../data/data-in/df_redlist.RDS")

# list of modelled species
sp_n_all <- unique(df_probs$sp_n)

# list of input species
colnames(df_data_in)[23:658]


# number of utmgrids with plant data
nbutm <- length(unique(df_data_in$utmID))

# probabilitie for Current scenario
df_probsCur <- filter(df_probs, projname == "allsp_Current")
remove(df_probs)

# Calculate prevalences
df_rare <- df_data_in %>%
  dplyr::select(sp_n_all) %>%
  gather(key = sp_n, value = presence) %>%
  group_by(sp_n) %>%
  summarize(rareness = 1 - (sum(presence, na.rm = TRUE) / nbutm))

df_rare <- df_rare %>%
  arrange(rareness) %>%
  mutate(ranksp = row_number() /  nrow(df_rare))

df_rare <- df_rare %>%
  mutate(top10 = ifelse(ranksp >= 0.9, 1, 0))

df_rare$exp <- 10000^(df_rare$rareness)
df_rare$exp2 <- df_rare$exp / 10000

ggplot(df_rare, aes(x = rareness, y = exp)) + geom_point()
ggplot(df_rare, aes(x = rareness, y = exp2)) + geom_point()
ggplot(df_rare, aes(x = rareness, y = ranksp)) + geom_point() +
  geom_point(aes(y = exp2, colour = "red"))

# Rode listsoorten weging
head(df_redlist)

temp <- filter(df_redlist, !is.na(RLweging))

temp2 <- inner_join(x = df_rare,
                  y = df_redlist %>%
                         dplyr::select(Soortcode, RLweging),
                  by = c("sp_n" = "Soortcode"))


temp3 <- filter(temp2, !is.na(RLweging))

# join prevalences to df_probs
df_probsCur2 <- df_probsCur %>%
  left_join(df_rare,
            by = "sp_n")

df_cur_sump <- df_probsCur2 %>%
  mutate(probsr = probs * top10,
         binsr = bins * top10) %>%
  group_by(utmID) %>%
  summarise(probsum = sum(probs) / 1000,
            binsum = sum(bins),
            probsumr = sum(probsr) / 1000,
            binsumr = sum(binsr))

df_cur_sump <- df_cur_sump %>%
  left_join(df_proj_in %>%
              dplyr::select(utmID, X, Y),
            by = "utmID")

df_cur_sump$probv <- df_cur_sump$probsum - df_cur_sump$probsumr

ggplot(df_cur_sump, aes(x = X, y = Y, colour = probsum)) +
  scale_colour_gradientn(colours = rainbow(3)) +
  geom_point() + coord_fixed()

ggplot(df_cur_sump, aes(x = X, y = Y, colour = probsumr)) +
  scale_colour_gradientn(colours = rainbow(3)) +
  geom_point() + coord_fixed()

ggplot(df_cur_sump, aes(x = X, y = Y, colour = probsum - probsumr)) +
  scale_colour_gradientn(colours = rainbow(3)) +
  geom_point() + coord_fixed()


ggplot(df_cur_sump, aes(x = probsum, y = probsumr)) + geom_point()

ggplot(df_cur_sump, aes(x = binsum, y = binsumr)) + geom_point()

ggplot(df_cur_sump, aes(x = X, y = Y, colour = binsumr)) +
  scale_colour_gradientn(colours = rainbow(3)) +
  geom_point() + coord_fixed()




write.table(df_cur_sump, file = paste0(path, "df_cur_sump.txt"),
            quote = FALSE, sep = ";", dec = ",", row.names = FALSE)

