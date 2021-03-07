
# load libraries

library("ggplot2")
library("scales")
library("grid")
library("vegan")
library("Biostrings")

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")

# load paths 

source("paths.R")

# files

design.file <- paste(data.dir, "design_IPL.txt", sep="")
otu_table.file <- paste(results.dir, "otu_table_IPL.txt", sep="")
rep_seqs.file <- paste(results.dir, "rep_seqs_IPL.fasta", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)
rep_seqs <- readDNAStringSet(rep_seqs.file)

# re-order data matrices

idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# normalize otu tables

otu_table_norm <- apply(otu_table, 2, function(x) x/sum(x))

# split OTU tables

idx <- design$system %in% c("greenhouse") & design$compartment %in% c("root")
otu_table_At_NC <- otu_table[, idx]
otu_table_norm_At_NC <- otu_table_norm[, idx]

idx <- design$system %in% c("IPL")
otu_table_CrIPL <- otu_table[, idx]
otu_table_norm_CrIPL <- otu_table_norm[, idx]

idx <- design$system %in% c("greenhouse") & design$compartment %in% c("phycosphere")
otu_table_Cr_NC <- otu_table[, idx]
otu_table_norm_Cr_NC <- otu_table_norm[, idx]

idx <- design$system %in% c("flask") & design$media %in% c("TP", "BD")
otu_table_Cr_MC <- otu_table[, idx]
otu_table_norm_Cr_MC <- otu_table_norm[, idx]

