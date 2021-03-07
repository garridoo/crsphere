
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

### Cooloola site dataset (Yeoh et. al., 2017)

# load paths 

source("paths.R")

# files

design.file <- paste(data.dir, "design_yeoh.txt", sep="")
asv_table.file <- paste(results.dir, "asv_table_yeoh.txt", sep="")
taxonomy.file <- paste(results.dir, "taxonomy_yeoh.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
asv_table <- read.table(asv_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, fill=T, quote="")

# re-order data matrices

idx <- design$SampleID %in% colnames(asv_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(asv_table))
asv_table <- asv_table[, idx]

idx <- match(rownames(asv_table), taxonomy[, 1])
taxonomy <- taxonomy[idx, ]

# keep only bacterial AVSs

idx <- taxonomy$kingdom=="Bacteria" & taxonomy$order!="Chloroplast"
taxonomy <- taxonomy[idx, ]
asv_table <- asv_table[idx, ]

# normalize otu tables

design$depth <- colSums(asv_table)
asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

# thresholding

threshold <- 0
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
# idx <- apply(asv_table_norm, 1, mean) * 100 >= threshold
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx  <- design$depth >= 500
design  <- design[idx, ]
asv_table  <- asv_table[, idx]
asv_table_norm  <- asv_table_norm[, idx]

# re-order taxonomy table

idx <- match(rownames(asv_table), taxonomy[, 1]) 
taxonomy <- taxonomy[idx, ]

