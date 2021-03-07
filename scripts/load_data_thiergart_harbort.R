
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# load paths 

source("paths.R")

# files

design.file <- paste(data.dir, "design_thiergart_harbort.txt", sep="")
asv_table.file <- paste(results.dir, "asv_table_thiergart_harbort.txt", sep="")
taxonomy.file <- paste(results.dir, "taxonomy_thiergart_harbort.txt", sep="")

# load data

design <- read.table(design.file, header=T, sep="\t")
asv_table <- read.table(asv_table.file, sep="\t", header=T, check.names=F)
taxonomy <- read.table(taxonomy.file, sep="\t", header=T, check.names=F)

# re-order data matrices

idx <- design$Original.SampleID %in% colnames(asv_table)
design <- design[idx, ]

idx <- match(design$Original.SampleID, colnames(asv_table))
asv_table <- asv_table[, idx]

colnames(asv_table) <- design$SampleID

# remove chloroplast reads

chloro <- taxonomy$ASV[taxonomy$Phylum=="Cyanobacteria/Chloroplast"]
asv_table <- asv_table[!rownames(asv_table) %in% chloro, ]

# normalize asv table

design$depth <- colSums(asv_table)
asv_table_norm <- apply(asv_table, 2, function(x) x/sum(x))

# thresholding

threshold <- 0
idx <- rowSums(asv_table_norm * 100 > threshold) >= 1
asv_table <- asv_table[idx, ]
asv_table_norm <- asv_table_norm[idx, ]

idx  <- design$depth >= 500
design  <- design[idx, ]
asv_table  <- asv_table[, idx]
asv_table_norm  <- asv_table_norm[, idx]

# re-order taxonomy table

idx <- match(rownames(asv_table), taxonomy[, 1]) 
taxonomy <- taxonomy[idx, ]

