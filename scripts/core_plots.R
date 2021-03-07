
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# load libraries

library(ape)

# calculate core taxonomic groups and their relative
# abundances across multiple host species and soil types

source("core.R")

# load paths 

source("paths.R")

# read plant rubisco (rbcL) tree (adapted and expanded from Yeoh et al., 2017)

tree.file <- paste(results.dir, "rbcL.newick", sep="")
tree <- read.tree(tree.file)
tree$tip.label <- gsub("_", " ", tree$tip.label)

# draw dotplot of mean relative abundances

df <- melt(core_table_norm)
colnames(df) <- c("core_taxon", "host_species", "RA")

# reorder host species according to the rbcL phylogeny

l <- tree$tip.label
l <- c(l, levels(df$host_species)[!levels(df$host_species) %in% l])
df$host_species <- factor(df$host_species, levels=l)

p1 <- ggplot(df, aes(x=core_taxon, y=host_species)) +
             geom_point(aes(size=RA), shape=21, color="grey", fill="grey") +
             labs(x="", y="") +
             main_theme +
             theme(axis.text.x=element_text(size=10)) +
             theme(axis.text.y=element_text(size=10)) +
             theme(axis.title=element_text(size=10)) +
             theme(axis.text.x=element_text(angle=45, hjust=1)) +
             theme(legend.position="left")
 
ggsave(paste(figures.dir, "core_dotplot.pdf", sep=""), p1, width=8, height=10)

p1 <- ggplot(df, aes(x=RA*100+1, y=host_species, fill=core_taxon)) +
             geom_bar(stat="identity", position=position_dodge()) +
             labs(x="", y="") +
             scale_x_log10(breaks=c(1, 2, 11, 51, 101), limits=c(NA, 101)) +
             main_theme +
             theme(axis.text.y=element_text(size=14, face="italic", color="black")) +
             theme(axis.text.x=element_text(size=14, color="black", angle=45, hjust=1)) +
             theme(legend.position="left")
 
ggsave(paste(figures.dir, "core_barplot.pdf", sep=""), p1, width=8, height=10)
    
