# Proportion stats for paper
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)

source("../set_paths.R")
setPaths()

proportions <- read.csv(paste0(other.tables.path, "processing_input/composition.csv"))

mydb <- dbConnect(SQLite(), paste0(other.tables.path, "HI_tables.sqlite"))
meta <- dbReadTable(mydb, "proc_metadata")
dbDisconnect(mydb)


# t-test purity, diabetes status
res.t2d <- t.test(puritypercentage ~ diagnosis, data = meta[meta$diagnosis %in% c("Type2", "None"),], var.equal = TRUE)
res.t1d <- t.test(puritypercentage ~ diagnosis, data = meta[meta$diagnosis %in% c("Type1", "None"),], var.equal = TRUE)

ggplot(meta, aes(x = puritypercentage, group = diagnosis, fill = diagnosis)) +
  geom_density(alpha = 0.2) +
  theme_bw()

# correlation between purity and non-endocrine proportion
cor.test(meta$puritypercentage, meta$exo_per, use = "complete")
