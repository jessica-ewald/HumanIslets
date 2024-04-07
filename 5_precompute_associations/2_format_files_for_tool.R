# create tool files for feature search
# Jessica Ewald
# November 8, 2023

library(RSQLite)
library(dplyr)

proc.summary <- read.csv("/Users/jessicaewald/Desktop/RestTest/resources/humanislets/proc_variable_summary.csv")

mydb <- dbConnect(SQLite(), "/Users/jessicaewald/sqlite/HI_precomputed.sqlite")
dbListTables(mydb)
omics_outcomes <- dbReadTable(mydb, "omics_outcomes")
dbDisconnect(mydb)


# create variable search file
genes <- omics_outcomes[,c("Feature", "Gene_ID", "Description")]
genes$Description <- gsub('\\"', '', genes$Description)
genes <- distinct(genes)
gene.entries <- data.frame(label = paste0(genes$Feature, " (Entrez: ", genes$Gene_ID, " - ", genes$Description, ")"),
                           type = "omics", ID = genes$Gene_ID, Contrast = NA)



outcomes <- omics_outcomes[,c("Metadata_ID", "Contrast", "Metadata_group")]
outcomes <- distinct(outcomes)

# pull out ephys outcomes
ephys.outcomes <- outcomes[outcomes$Metadata_group == "Single-cell Function", ]
outcomes <- outcomes[outcomes$Metadata_group != "Single-cell Function", ]

# add other info to normal variables
outcomes <- merge(outcomes, proc.summary, by.x = "Metadata_ID", by.y = "column", all.x = TRUE, all.y = FALSE)

# pull out outcomes with a contrast
outcomes.contrast <- outcomes[!is.na(outcomes$Contrast),]
outcomes <- outcomes[is.na(outcomes$Contrast), ]

# process ephys outcomes
ephys.outcomes <- ephys.outcomes[,c("Metadata_ID", "Metadata_group")]
ephys.outcomes$second <- gsub(".*_cell", "", ephys.outcomes$Metadata_ID)
ephys.outcomes$second <- gsub(".*_donor", "", ephys.outcomes$second)

ephys.outcomes$column <- ephys.outcomes$Metadata_ID
ephys.outcomes$Metadata_ID <- gsub("_cell.*", "", ephys.outcomes$Metadata_ID)
ephys.outcomes$Metadata_ID <- gsub("_donor.*", "", ephys.outcomes$Metadata_ID)

ephys.outcomes$Metadata_group[grep("_cell_", ephys.outcomes$column)] <- "Cell-level"
ephys.outcomes$Metadata_group[grep("_donor_", ephys.outcomes$column)] <- "Donor-level"
ephys.outcomes$second <- sub('.', '', ephys.outcomes$second)
ephys.outcomes$cell <- gsub("_.*", "", ephys.outcomes$second)
ephys.outcomes$gluc <- gsub(".*_", "", ephys.outcomes$second)

ephys.names <- proc.summary[proc.summary$group_id == "ephys_cell", c("column", "display")]
ephys.names$column <- gsub("_cell", "", ephys.names$column)
ephys.outcomes <- merge(ephys.outcomes, ephys.names, by.x = "Metadata_ID", by.y = "column")

# format outcome entries
outcome.names <- data.frame(label = paste0(outcomes$display, " [", outcomes$group, "]"), 
                            type = "outcome", ID = outcomes$Metadata_ID, Contrast = NA)

outcome.names.contrast <- data.frame(label = paste0(outcomes.contrast$display, " [", outcomes.contrast$Contrast, ", ", outcomes.contrast$group, "]"),
                                     type = "outcome", ID = outcomes.contrast$Metadata_ID, Contrast = outcomes.contrast$Contrast)

ephys.outcome.names <- data.frame(label = paste0(ephys.outcomes$display, " [", ephys.outcomes$cell, " cells & ", ephys.outcomes$gluc, " glucose, ",
                                                 ephys.outcomes$Metadata_group, " Single-cell Function]"),
                                  type = "outcome_ephys", ID = ephys.outcomes$column, Contrast = NA)
  
# combine all entries into one data.frame
all.names <- rbind(outcome.names, outcome.names.contrast)
all.names <- rbind(all.names, ephys.outcome.names)
all.names <- rbind(all.names, gene.entries)

# write out results
write.csv(all.names, "/Users/jessicaewald/Desktop/RestTest/resources/humanislets/feature_search_names.csv", row.names = FALSE)
