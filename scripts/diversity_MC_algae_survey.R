
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load At Cr natural community dataset

source("load_data_MC_algae_survey.R")

# subset samples of interest

idx <- (design$timepoint %in% c("14", "21", "35") | (design$timepoint=="0" & design$community %in% c("Input", "Unplanted"))) & 
       design$community %in% c("Unplanted", "Chlamydomonas reinhardtii", "Chlamydomonas sp. (1030)", "Spirogloea muscicola", "Microthamnion sp. (1108)", "Klebsormidium sp. (1121)") &
       design$environment %in% c("Input", "Unplanted", "terrestrial", "aquatic") &
       design$medium %in% c("SFW") &
       T
    
design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]

### beta diversity

colors <- data.frame(group=c("Chlamydomonas reinhardtii", "Chlamydomonas sp. (1030)", "Spirogloea muscicola", "Microthamnion sp. (1108)", "Klebsormidium sp. (1121)", "Unplanted"),
                     color=c(cr_color,                    chlamy_1030_color,          spirogloea_color,       microthamnion_color,         klebsormidium_1121_color,   soil_color))

shapes <- data.frame(group=c("Unplanted", "terrestrial",     "aquatic"),
                     shape=c(soil_shape,  phycosphere_shape, aquatic_shape))

### CPCoA analysis

colors_cpcoa <- data.frame(group=c("Chlamydomonas reinhardtii", "Chlamydomonas sp. (1030)", "Spirogloea muscicola", "Microthamnion sp. (1108)", "Klebsormidium sp. (1121)", "Unplanted"),
                           color=c(cr_color,                    chlamy_1030_color,          spirogloea_color,       microthamnion_color,         klebsormidium_1121_color,   soil_color))

shapes_cpcoa <- data.frame(group=c("Unplanted", "terrestrial",     "aquatic"),
                           shape=c(soil_shape,  phycosphere_shape, aquatic_shape))

sqrt_transform <- T

capscale.gen <- capscale(bray_curtis ~ community * environment + Condition(replicate),
                         data=design_subset, add=F, sqrt.dist=sqrt_transform)

# ANOVA-like permutation analysis

perm_anova.gen <- anova.cca(capscale.gen)
print(perm_anova.gen)
                                                    
# generate variability tables and calculate confidence intervals for the variance

var_tbl.gen <- variability_table(capscale.gen)

eig <- capscale.gen$CCA$eig

variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi

variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[, 1:2]
colnames(points_cpcoa) <- c("x", "y")

points_cpcoa <- cbind(points_cpcoa, design_subset[match(rownames(points_cpcoa), design_subset$SampleID), ])

colors_cpcoa <- colors_cpcoa[colors_cpcoa$group %in% points_cpcoa$community, ]
points_cpcoa$community <- factor(points_cpcoa$community, levels=colors_cpcoa$group)

shapes_cpcoa <- shapes_cpcoa[shapes_cpcoa$group %in% points_cpcoa$environment, ]
points_cpcoa$environment <- factor(points_cpcoa$environment, levels=shapes_cpcoa$group)

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ community, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('community','seg_x','seg_y')), by='community', sort=FALSE)

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=community, shape=environment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none") # + theme(aspect.ratio=1)

ggsave(paste(figures.dir, "algae_survey_MC_subset_CPCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### Hierarchical clustering dendrogram and heatmap

df <- melt(asv_table_subset)
colnames(df) <- c("asv", "SampleID", "ra")
idx <- match(df$SampleID, design_subset$SampleID)
df <- cbind(df, design_subset[idx, -1])

threshold <- 0.1
idx <- apply(asv_table_subset, 1, FUN=mean) * 100 > threshold
m <- asv_table_subset[idx, ]

d <- vegdist(t(asv_table_subset), method="bray")
cluster = hclust(d, method="ward.D2")

m <- log(m*1000+1)

colors_hc <- data.frame(group=c("Chlamydomonas reinhardtii", "Chlamydomonas sp. (1030)", "Spirogloea muscicola", "Microthamnion sp. (1108)", "Klebsormidium sp. (1121)", "Unplanted"),
                        color=c(cr_color,                    chlamy_1030_color,          spirogloea_color,       microthamnion_color,         klebsormidium_1121_color,   soil_color))

idx <- match(design_subset$community, colors_hc$group)
colcols <- as.character(colors_hc$color[idx])

pdf(paste(figures.dir, "algae_survey_NC_subset_heatmap.pdf", sep=""), p, width=10, height=10)
p <- heatmap.2(m,
               Colv=as.dendrogram(cluster),
               ColSideColors=colcols,
               trace="none",
               labRow=FALSE, labCol=FALSE)
dev.off()

