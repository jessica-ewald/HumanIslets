# Process Raw RNA-seq counts for web-tool
# Jessica Ewald
# April 25, 2023

library(dplyr)
library(ggplot2)
library(RSQLite)
library(edgeR)
library(sva)

# read in 
fileName <- "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Omics_data/raw/raw_rnaseq.txt"
counts <- data.table::fread(fileName, header=TRUE, check.names=FALSE, data.table=FALSE)

# annotate to Entrez
myDb <- dbConnect(SQLite(), "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/annotation_libraries/hsa_genes.sqlite")
ensg <- dbReadTable(myDb, "entrez_embl_gene")
entrez <- dbReadTable(myDb, "entrez")
dbDisconnect(myDb)

feature.vec <- counts$V1
hit.inx <- match(feature.vec, ensg[, "accession"])
gene.ids <- ensg[hit.inx, ]
counts$V1 <- gene.ids$gene_id
counts <- counts[!is.na(counts$V1), ]

# remove duplicates
counts <- aggregate(counts[-1], counts[1], sum)
rownames(counts) <- counts$V1
counts <- counts[,-1]

# remove 80% zeros
num.zeros <- apply(counts, 1, function(x){sum(x == 0)})
inds.remove <- which(num.zeros > 0.8*dim(counts)[2])
counts <- counts[-inds.remove, ]

# normalize for sequencing depth
nf <- calcNormFactors(counts, method="RLE")
lcpm <- cpm(counts, lib.size=colSums(counts)*nf, log = TRUE)

# PCA before batch effect
batch <- read.table("/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Omics_data/orig/HumanIslets_OxfordStanfordCombined_GeneExp_SampleInfo.txt",
                    sep = "\t")
colnames(batch) <- c("record_id", "batch")

pr <- prcomp(t(lcpm), scale = TRUE)
pr.x <- pr$x
pr.x <- merge(pr.x, batch, by.x = "row.names", by.y = "record_id")

ggplot(pr.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("LCPM")

# one batch (Batch2_Oxford) is EXTREMELY different. remove and re-normalize
batch <- batch[batch$batch != "Batch2_Oxford",]
counts <- counts[,colnames(counts) %in% batch$record_id]
nf <- calcNormFactors(counts, method="RLE")
lcpm <- cpm(counts, lib.size=colSums(counts)*nf, log = TRUE)

pr <- prcomp(t(lcpm), scale = TRUE)
pr.x <- pr$x
pr.x <- merge(pr.x, batch, by.x = "row.names", by.y = "record_id", all = F)

ggplot(pr.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("LCPM (Oxford2 batch removed)")

# adjust for batch effect
identical(batch$record_id, colnames(counts)) # check to make sure order preserved
counts <- as.matrix(counts)
adjusted <- ComBat_seq(counts, batch=batch$batch, group=NULL)

nf <- calcNormFactors(adjusted, method="RLE")
lcpm <- cpm(adjusted, lib.size=colSums(adjusted)*nf, log = TRUE)

pr <- prcomp(t(lcpm), scale = TRUE)
pr.x <- pr$x
pr.x <- merge(pr.x, batch, by.x = "row.names", by.y = "record_id", all = F)

ggplot(pr.x, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  xlab(paste0("PC1 (", round((pr$sdev[1]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  ylab(paste0("PC2 (", round((pr$sdev[2]^2)/sum(pr$sdev^2)*100,1), "%)")) +
  theme_bw() +
  ggtitle("LCPM (adjusted for batch effect)")

lcpm[lcpm == min(lcpm)] <- NA
lcpm[lcpm == min(lcpm, na.rm = TRUE)] <- NA

# filter based on variance
lcpm.var <- apply(lcpm, 1, var, na.rm = TRUE)
lcpm.mean <- apply(lcpm, 1, mean, na.rm = TRUE)
plot(lcpm.mean, lcpm.var)
var.thresh <- quantile(lcpm.var, 0.20)
lcpm.keep <- which(lcpm.var > var.thresh)
lcpm <- lcpm[lcpm.keep, ]

# # only keep samples with >30% purity
# mydb <- dbConnect(SQLite(), "/Users/jessicaewald/sqlite/HI_tables.sqlite")
# iso <- dbReadTable(mydb, "isolation")
# dbDisconnect(mydb)
# 
# high.purity <- iso$record_id[iso$puritypercentage > 30]
# lcpm <- lcpm[,colnames(lcpm) %in% high.purity]

# write out files
write.csv(adjusted, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Omics_data/raw/raw_rnaseq_batch.csv")
write.csv(lcpm, "/Users/jessicaewald/Library/CloudStorage/OneDrive-McGillUniversity/XiaLab/Tools/HumanIslets/Omics_data/proc/proc_rnaseq.csv")


