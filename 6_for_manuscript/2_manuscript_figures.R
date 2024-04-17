# Figures for manuscript
# Author: Jessica Ewald

## Set your working directory to the "6_for_manuscript" directory

library(dplyr)
library(RSQLite)

source("../set_paths.R")
setPaths()


###### Figure 2A #######
# This figure shows the distribution of the main metadata across all donors

# get data
mydb <- dbConnect(RSQLite::SQLite(), paste0(sqlite.path, "HI_tables.sqlite"))
donor <- dbReadTable(mydb, "donor")
dbDisconnect(mydb)

# sex
pdf(file="./figures/fig2A-1.pdf", width=4, heigh=4)
ggplot(donor, aes(x = donorsex)) +
  geom_bar(stat="count") +
  theme_classic() +
  xlab("Sex") +
  ylab("Count")
dev.off()

# diagnosis
pdf(file="./figures/fig2A-2.pdf", width=4, heigh=4)
ggplot(donor, aes(x = diagnosis)) +
  geom_bar(stat="count") +
  theme_classic() +
  xlab("Diabetes diagnosis") +
  ylab("Count")
dev.off()