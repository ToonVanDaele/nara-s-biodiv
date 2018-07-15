## Plot explanatory variables

library(tidyverse)

modelname <- "rlw"
path <- paste0("../data/data-out/", modelname, "/")

df_expl <- readRDS("../data/data-in/df_expl_eda.RDS")
explspec <- readRDS(paste0(path, "explspec.RDS"))

vars <- unique(explspec$varname)

df_expl_spr <- df_expl %>%
  filter(!varname == "Opp_m2") %>%
  split(list(.$scen, .$varname)) %>%
  walk(function(df) {
        title <- paste0(unique(df$varname), "_", unique(df$scen))
        p <- ggplot(df, aes(x = X, y = Y, colour = value)) +
          geom_point(size = 0.6) + coord_fixed() + ggtitle(title)
        ggsave(filename = paste0(path, "z_", title, ".png"), p)})

