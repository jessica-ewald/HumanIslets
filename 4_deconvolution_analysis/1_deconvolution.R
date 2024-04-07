# Deconvolution version 2
# Jessica Ewald
# Nov 30, 2023

library(dplyr)
library(stringr)
library(RSQLite)
library(ggplot2)
library(pheatmap)
library(ggpubr)

# read in proteomics data
prot <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Omics_data/raw/raw_prot.csv",
                 header = TRUE, row.names = 1)

# process uniprot annotation data
uniprot.dat <- read.table("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Rscripts/deconvolution/proteomics/uniprot_id_mass.tsv",
                          header = TRUE, sep = "\t", fill = TRUE, quote="")
IDs <- str_split(uniprot.dat$Gene.Names, " ")
names(IDs) <- uniprot.dat$Entry
IDs <- data.table::melt(IDs)

# normalize by sample volume
sample.mass <- apply(prot, 2, sum, na.rm = TRUE)
prot <- sweep(prot, 2, sample.mass, "/")
prot[is.na(prot)] <- min(prot, na.rm = T)/5
  
########## COMPUTE CELL TYPE PROPORTIONS ###########

# purity & cell type prior
max.purity <- 0.98
proB.mean <- 0.55
proA.mean <- 0.30
proD.mean <- 0.10
proG.mean <- 0.05

# Read in and process marker genes
markers <- readRDS("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Rscripts/deconvolution/proteomics/markers_herrera.rds")

# get marker genes
bcm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "beta"], ]
acm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "alpha"], ]
dcm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "delta"], ]
gcm.dat <- prot[rownames(prot) %in% markers$uniprot[markers$cell_type == "gamma"], ]

# compute median of each marker
bcell <- apply(bcm.dat, 1, median, na.rm = TRUE) %>% unlist()
bcell <- data.frame(uniprot = names(bcell), median.bulk = unname(bcell))

acell <- apply(acm.dat, 1, median, na.rm = TRUE) %>% unlist()
acell <- data.frame(uniprot = names(acell), median.bulk = unname(acell))

dcell <- apply(dcm.dat, 1, median, na.rm = TRUE) %>% unlist()
dcell <- data.frame(uniprot = names(dcell), median.bulk = unname(dcell))

gcell <- apply(gcm.dat, 1, median, na.rm = TRUE) %>% unlist()
gcell <- data.frame(uniprot = names(gcell), median.bulk = unname(gcell))


# compute estimated median level of marker within cell type
bcell$median.beta <- bcell$median.bulk/proB.mean
acell$median.alpha <- acell$median.bulk/proA.mean
dcell$median.delta <- dcell$median.bulk/proD.mean
gcell$median.gamma <- gcell$median.bulk/proG.mean

# re-order
bcm.dat <- bcm.dat[match(bcell$uniprot, rownames(bcm.dat)), ]
acm.dat <- acm.dat[match(acell$uniprot, rownames(acm.dat)), ]
dcm.dat <- dcm.dat[match(dcell$uniprot, rownames(dcm.dat)), ]
gcm.dat <- gcm.dat[match(gcell$uniprot, rownames(gcm.dat)), ]

# check
identical(bcell$uniprot, rownames(bcm.dat)) & 
  identical(acell$uniprot, rownames(acm.dat)) & 
  identical(dcell$uniprot, rownames(dcm.dat)) & 
  identical(gcell$uniprot, rownames(gcm.dat))

# compute proportion of each cell type in each sample for each marker gene
bcell.pro <- sweep(bcm.dat, 1, bcell$median.beta, "/")
acell.pro <- sweep(acm.dat, 1, acell$median.alpha, "/")
dcell.pro <- sweep(dcm.dat, 1, dcell$median.delta, "/")
gcell.pro <- sweep(gcm.dat, 1, gcell$median.gamma, "/")

# compute median proportion of each cell type in each sample, across all marker gene estimates
bcell.pro.med <- apply(bcell.pro, 2, median, na.rm = TRUE)
acell.pro.med <- apply(acell.pro, 2, median, na.rm = TRUE)
dcell.pro.med <- apply(dcell.pro, 2, median, na.rm = TRUE)
gcell.pro.med <- apply(gcell.pro, 2, median, na.rm = TRUE)

# check
identical(names(bcell.pro.med), names(acell.pro.med)) &
  identical(names(bcell.pro.med), names(dcell.pro.med)) &
  identical(names(bcell.pro.med), names(gcell.pro.med))

# compute proportions
cell.pro <- data.frame(record_id = names(acell.pro.med), beta = bcell.pro.med, alpha = acell.pro.med, 
                       delta = dcell.pro.med, gamma = gcell.pro.med)
cell.pro$end_total <- apply(cell.pro[,-1], 1, sum)

# Assume max possible purity is 98%. Calculate the final "total" from which to subtract everything else.
max.tot.end <- max(cell.pro$end_total)
mass.total <- max.tot.end/max.purity
cell.pro$exo_total <- mass.total - cell.pro$end_total
cell.pro$total <- cell.pro$end_total + cell.pro$exo_total

cell.pro$beta_per <- cell.pro$beta/cell.pro$total
cell.pro$alpha_per <- cell.pro$alpha/cell.pro$total
cell.pro$delta_per <- cell.pro$delta/cell.pro$total
cell.pro$gamma_per <- cell.pro$gamma/cell.pro$total
cell.pro$exo_per <- cell.pro$exo_total/cell.pro$total

cell.pro$beta_end <- cell.pro$beta/cell.pro$end_total
cell.pro$alpha_end <- cell.pro$alpha/cell.pro$end_total
cell.pro$delta_end <- cell.pro$delta/cell.pro$end_total
cell.pro$gamma_end <- cell.pro$gamma/cell.pro$end_total

saveRDS(cell.pro, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Rscripts/deconvolution/proteomics/v3_results/proportions.rds")
write.csv(cell.pro[,c("record_id", "beta_end", "alpha_end", "delta_end", "gamma_end", "exo_per")],
          "/Users/jessicaewald/Desktop/RestTest/resources/humanislets/processing_input/composition.csv",
          row.names = FALSE)

########## ANALYZE CELL TYPE PROPORTIONS ###########

library(GGally)
ggpairs(cell.pro[,c(9:13)], aes(alpha = 0.4)) + theme_bw() 
ggpairs(cell.pro[,c(14:17)], aes(alpha = 0.4)) + theme_bw() 

# plot proportions
stacked.bar <- cell.pro[,c(1:5,7)]

stacked.bar <- reshape2::melt(stacked.bar)
colnames(stacked.bar) <- c("record_id", "cell_type", "proportion")
stacked.bar$cell_type <- as.character(stacked.bar$cell_type)
stacked.bar$cell_type[stacked.bar$cell_type == "exo_total"] <- "non.endocrine"
stacked.bar$record_id <- factor(stacked.bar$record_id, levels = cell.pro$record_id[order(cell.pro$exo_total, cell.pro$beta)])

ggplot(stacked.bar, aes(fill=cell_type, y=proportion, x=record_id)) + 
  geom_bar(position='stack', stat='identity', width=1) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6, margin = margin(r = -15)),
        panel.grid.major.y = element_blank()) +
  coord_flip()

# correlation structure between different cell type measures
cormat <- cor(cell.pro[,-1])
cormat[lower.tri(cormat)] <- NA
cormat <- round(cormat, 2) %>% reshape2::melt(., na.rm = TRUE)


# Plot proportions
# Like Fig3A in this paper: https://www.nature.com/articles/s41598-017-16300-w 

prop <- cell.pro[,c("record_id", "beta_end", "alpha_end", "delta_end", "gamma_end")]
prop <- prop[order(prop$beta_end, decreasing = TRUE), ]
prop$record_id <- factor(prop$record_id, levels = prop$record_id)
prop <- reshape2::melt(prop)
colnames(prop) <- c("record_id", "Cell_Type", "Proportion")

ggplot(prop, aes(x = record_id, y = Proportion, color = Cell_Type, group = Cell_Type)) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))

# measured proportions
meas <- read.csv("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Rscripts/deconvolution/cell_proportions_HIPP.csv")
colnames(meas) <- c("record_id", "beta.m", "alpha.m", "delta.m", "exo.m", "endo.m")

prop <- cell.pro[,c("record_id", "alpha", "beta", "delta")]
prop$alpha <- prop$alph
prop$beta <- prop$beta
prop$delta <- prop$delta

prop$m.total <- apply(prop[,2:4], 1, sum)
prop$beta.md <- prop$beta/prop$m.total*100
prop$alpha.md <- prop$alpha/prop$m.total*100
prop$delta.md <- prop$delta/prop$m.total*100

compare <- merge(prop[,c("record_id", "beta.md", "alpha.md", "delta.md")], meas, by = "record_id")

# there is a correlation and all p-values are significant
# could it agree better if a surface area / volume correction were applied?
# proportion volume/mass, proportion cell number, proportion stained pixels (sort of surface area??)
# given systematically different cell sizes & shapes between types, these would all give different 
# numbers

p <- ggplot(compare, aes(x = beta.m, y = beta.md)) +
  geom_point() +
  theme_bw()
p + geom_abline(intercept = 0, slope = 1)

cor.test(compare$beta.m, compare$beta.md)

p <- ggplot(compare, aes(x = alpha.m, y = alpha.md)) +
  geom_point() +
  theme_bw()
p + geom_abline(intercept = 0, slope = 1)
cor.test(compare$alpha.m, compare$alpha.md)

p <- ggplot(compare, aes(x = delta.m, y = delta.md)) +
  geom_point() +
  theme_bw()
p + geom_abline(intercept = 0, slope = 1)
cor.test(compare$delta.m, compare$delta.md)




