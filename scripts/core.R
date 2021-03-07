
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

options(warn=-1)

# cleanup

rm(list=ls())

# load plotting functions

source("plotting_functions.R")
source("plotting_parameters.R")
source("cpcoa.func.R")

# load plotting functions

library("ggplot2")
library("scales")
library("grid")
library("vegan")

# ancillary  functions

core <- function(x) {
    
    # returns percentage of samples in which an feature is found

    if (sum(x)==0) return(0)
    else return(as.numeric(sum(x>0)/length(x)))

}

# occupancy threshold (core definition)

occupancy_threshold <- 0.6

### Cooloola site dataset (Yeoh et al., 2017)

source("load_data_yeoh.R")

# aggregate relative abundances to the order level

order_table_norm <- aggregate(asv_table_norm, by=list(taxonomy$order), FUN=sum)
rownames(order_table_norm) <- order_table_norm[, 1]
order_table_norm <- order_table_norm[, -1]

# subset samples of interest

idx <- design$type %in% c("root") &
       TRUE

design_yeoh <- design[idx, ]
asv_table_yeoh <- asv_table[, idx]
asv_table_norm_yeoh <- asv_table_norm [, idx]
order_table_norm_yeoh <- order_table_norm[, idx]
taxonomy_yeoh <- taxonomy

# aggregate (mean) order RA table per host species 

os_table_norm_yeoh <- aggregate(t(order_table_norm_yeoh), by=list(design_yeoh$host_sp), FUN=mean)
rownames(os_table_norm_yeoh) <- os_table_norm_yeoh[, 1]
os_table_norm_yeoh <- t(os_table_norm_yeoh[, -1])

# calculate occupancy table

occupancy_table_yeoh <- aggregate(t(order_table_norm_yeoh), by=list(design_yeoh$host_sp), FUN=core)
rownames(occupancy_table_yeoh) <- occupancy_table_yeoh[, 1]
occupancy_table_yeoh <- t(occupancy_table_yeoh[, -1])

# calculate core orders (allow for one host spp. missing)

idx <- rowSums(occupancy_table_yeoh < occupancy_threshold)<=1
core_yeoh <- rownames(occupancy_table_yeoh)[idx]

# remove taxonomic groups without classification

idx <- !grepl("Incertae", core_yeoh) & core_yeoh!="V1"
core_yeoh <- core_yeoh[idx]

### At/Lj natural community datasets (Thiergart et al., 2019; Harbort et al., 2020)

source("load_data_thiergart_harbort.R")

# aggregate relative abundances to the order level

order_table_norm <- aggregate(asv_table_norm, by=list(taxonomy$order), FUN=sum)
rownames(order_table_norm) <- order_table_norm[, 1]
order_table_norm <- order_table_norm[, -1]

# subset samples of interest

idx <- design$Study %in% c("zgadzaj_2018", "masayoshi_2016") &
       design$Soil.Batch %in% c("CAS", "CAS10", "CAS11b") &
       design$Compartment %in% c("root") &
       TRUE

design_wippel <- design[idx, ]
asv_table_wippel <- asv_table[, idx]
asv_table_norm_wippel <- asv_table_norm[, idx]
order_table_norm_wippel <- order_table_norm[, idx]
taxonomy_wippel <- taxonomy

# aggregate (mean) order RA table per host species 

os_table_norm_wippel <- aggregate(t(order_table_norm_wippel), by=list(design_wippel$Host.Species), FUN=mean)
rownames(os_table_norm_wippel) <- os_table_norm_wippel[, 1]
os_table_norm_wippel <- t(os_table_norm_wippel[, -1])
colnames(os_table_norm_wippel) <- c("Arabidopsis thaliana (CAS 11)", "Lotus japonicus")

# calculate occupancy table

occupancy_table_wippel <- aggregate(t(order_table_norm_wippel), by=list(design_wippel$Host.Species), FUN=core)
rownames(occupancy_table_wippel) <- occupancy_table_wippel[, 1]
occupancy_table_wippel <- t(occupancy_table_wippel[, -1])

# calculate core orders

idx <- rowSums(occupancy_table_wippel < occupancy_threshold)==0
core_wippel <- rownames(occupancy_table_wippel)[idx]

# remove taxonomic groups with without classification

idx <- !grepl("Incertae", core_wippel) & core_wippel!="V1"
core_wippel <- core_wippel[idx]

### At/Cr natural community dataset (Duran et al., 2020)

source("load_data_NC.R")

# aggregate relative abundances to the order level

order_table_norm <- aggregate(asv_table_norm, by=list(taxonomy$order), FUN=sum)
rownames(order_table_norm) <- order_table_norm[, 1]
order_table_norm <- order_table_norm[, -1]

# subset samples of interest

idx <- design$system=="soil" &
       design$timepoint %in% c("Day36") &
       design$compartment %in% c("root", "phycosphere") &
       TRUE

design_NC <- design[idx, ]
asv_table_NC <- asv_table[, idx]
asv_table_norm_NC <- asv_table_norm[, idx]
order_table_norm_NC <- order_table_norm[, idx]
taxonomy_NC <- taxonomy

# aggregate (mean) order RA table per host species 

os_table_norm_NC <- aggregate(t(order_table_norm_NC), by=list(design_NC$compartment), FUN=mean)
rownames(os_table_norm_NC) <- os_table_norm_NC[, 1]
os_table_norm_NC <- t(os_table_norm_NC[, -1])
colnames(os_table_norm_NC) <- c("Chlamydomonas reinhardtii", "Arabidopsis thaliana")

# calculate occupancy table

occupancy_table_NC <- aggregate(t(order_table_norm_NC), by=list(design_NC$compartment), FUN=core)
rownames(occupancy_table_NC) <- occupancy_table_NC[, 1]
occupancy_table_NC <- t(occupancy_table_NC[, -1])

# calculate core orders

idx <- rowSums(occupancy_table_NC < occupancy_threshold)==0
core_NC <- rownames(occupancy_table_NC)[idx]

# remove taxonomic groups with without classification

idx <- !grepl("Incertae", core_NC) & core_NC!="V1"
core_NC <- core_NC[idx]

### Cr mesocosm dataset (Duran et al., 2020)

source("load_data_MC.R")

# aggregate relative abundances to the order level

order_table_norm <- aggregate(asv_table_norm, by=list(taxonomy$order), FUN=sum)
rownames(order_table_norm) <- order_table_norm[, 1]
order_table_norm <- order_table_norm[, -1]

# subset samples of interest

idx <- design$system=="flask" &
       design$community_AP %in% c("B+C_FALSE") &
       design$soil %in% c("CAS", "Golm") &
       design$experiment %in% c("AL1", "AL2", "AL3", "AL4", "AL5", "AL6") &
       design$medium %in% c("TP", "TP+AP") &
       design$timepoint %in% c("Day7") &
       TRUE

design_MC <- design[idx, ]
asv_table_MC <- asv_table[, idx]
asv_table_norm_MC <- asv_table_norm[, idx]
order_table_norm_MC <- order_table_norm[, idx]
taxonomy_MC <- taxonomy

# aggregate (mean) order RA table per host species 

os_table_norm_MC <- aggregate(t(order_table_norm_MC), by=list(design_MC$soil), FUN=mean)
rownames(os_table_norm_MC) <- os_table_norm_MC[, 1]
os_table_norm_MC <- t(os_table_norm_MC[, -1])
colnames(os_table_norm_MC) <- c("Chlamydomonnas reinhardtii (CAS)", "Chlamydomonnas reinhardtii (Golm)")

### embryophyte / chlorophyte core taxonomic groups

core_embryophyte <- intersect(core_yeoh, core_wippel)
core <- intersect(core_embryophyte, core_NC)

# get core taxa / host species average relative abundances tables

core_table_norm_yeoh <- os_table_norm_yeoh[match(core, rownames(os_table_norm_yeoh)), ]
core_table_norm_wippel <- os_table_norm_wippel[match(core, rownames(os_table_norm_wippel)), ]
core_table_norm_NC <- os_table_norm_NC[match(core, rownames(os_table_norm_NC)), ]
core_table_norm_MC <- os_table_norm_MC[match(core, rownames(os_table_norm_MC)), ]

core_table_norm <- cbind(core_table_norm_yeoh, core_table_norm_wippel, core_table_norm_NC, core_table_norm_MC)

