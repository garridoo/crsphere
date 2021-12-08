
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load At Cr natural community dataset

source("load_data_NC_nonflooded.R")

# subset samples of interest

idx <- (design$compartment %in% c("Input", "soil_surface", "phycosphere", "rhizosphere", "root") | (design$compartment=="soil" & design$timepoint==0)) &
       design$timepoint %in% c(0, 36, 47) &
       T
    
design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]

taxon_table <- aggregate(asv_table_subset, by=list(taxonomy$order), FUN=sum)
rownames(taxon_table) <- taxon_table[, 1]
taxon_table <- taxon_table[, -1]
asv_table_subset_b <- asv_table_subset

### beta diversity

design_subset$compartment <- as.character(design_subset$compartment)
idx <- design_subset$compartment %in% c("soil_surface") & design_subset$treatment %in% c("uncovered")
design_subset$compartment[idx] <- "CAS phycosphere"
idx <- design_subset$compartment %in% c("soil_surface") & design_subset$treatment %in% c("covered")
design_subset$compartment[idx] <- "soil"

colors <- data.frame(group=c("Input",    "soil",        "CAS phycosphere", "phycosphere", "rhizosphere", "root"),
                     color=c(soil_color, covered_color, uncovered_color,    cr_color,      at_color,      at_color))

shapes <- data.frame(group=c("Input",    "soil",        "CAS phycosphere", "phycosphere",     "rhizosphere",     "root"),
                     shape=c(soil_shape, soil_shape,    soil_shape,        phycosphere_shape, rhizosphere_shape, root_shape))

### CPCoA analysis

colors_cpcoa <- data.frame(group=c("Input",    "soil",        "CAS phycosphere", "phycosphere", "rhizosphere", "root"),
                     color=c(soil_color, covered_color, uncovered_color,    cr_color,      at_color,      at_color))

shapes_cpcoa <- data.frame(group=c("Input",    "soil",        "CAS phycosphere", "phycosphere",     "rhizosphere",     "root"),
                     shape=c(soil_shape, soil_surface_shape,    soil_surface_shape,        phycosphere_shape, rhizosphere_shape, root_shape))

design_subset$group <- as.character(design_subset$compartment)
design_subset$group[design_subset$compartment %in% c("phycosphere", "CAS phycosphere")] <- "phycosphere"


sqrt_transform <- T

capscale.gen <- capscale(bray_curtis ~ group,
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

points_cpcoa <- capscale.gen$CCA$wa[, 1:3]
colnames(points_cpcoa) <- c("x", "y", "z")

points_cpcoa <- cbind(points_cpcoa, design_subset[match(rownames(points_cpcoa), design_subset$SampleID), ])

colors_cpcoa <- colors_cpcoa[colors_cpcoa$group %in% points_cpcoa$compartment, ]
points_cpcoa$compartment <- factor(points_cpcoa$compartment, levels=colors_cpcoa$group)

shapes_cpcoa <- shapes_cpcoa[shapes_cpcoa$group %in% points_cpcoa$compartment, ]
points_cpcoa$compartment <- factor(points_cpcoa$compartment, levels=shapes_cpcoa$group)

# plot CPCo 1 and 2

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ compartment, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('compartment','seg_x','seg_y')), by='compartment', sort=FALSE)

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=compartment, shape=compartment, fill=compartment, label=SampleID)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_fill_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none") +
     theme(aspect.ratio=1)

ggsave(paste(figures.dir, "nonflooded_NC_CPCoA12.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

# plot CPCo 2 and 3

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$y, points_cpcoa$z) ~ compartment, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('compartment','seg_x','seg_y')), by='compartment', sort=FALSE)

p <- ggplot(points_cpcoa, aes(x=y, y=z, color=compartment, shape=compartment, fill=compartment, label=SampleID)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     # geom_text() +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_fill_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none") +
     theme(aspect.ratio=1)

ggsave(paste(figures.dir, "nonflooded_NC_CPCoA23.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)


# plot contribution of the core taxa to each compartment

# aggregate RAs of ASVs belonging to core taxa

core <- c("Caulobacterales", "Rhizobiales", "Sphingomonadales", "Burkholderiales", "Xanthomonadales", "Chitinophagales")
order_table_subset <- aggregate(asv_table_subset, by=list(taxonomy$order), FUN=sum)
rownames(order_table_subset) <- order_table_subset[, 1]
order_table_subset <- order_table_subset[, -1]

core_table_subset <- order_table_subset[rownames(order_table_subset) %in% core, ]
design_subset$core_RA <- colSums(core_table_subset)

# define shapes and colors

colors_core <- data.frame(group=c("Input",    "soil",        "CAS phycosphere", "phycosphere", "rhizosphere", "root"),
                          color=c(soil_color, covered_color, uncovered_color,   cr_color,      at_color,      at_color))

shapes_core <- data.frame(group=c("Input",    "soil",        "CAS phycosphere", "phycosphere",     "rhizosphere",     "root"),
                          shape=c(soil_shape, soil_shape,    soil_shape,        phycosphere_shape, rhizosphere_shape, root_shape))

df <- data.frame(compartment=design_subset$compartment, core_RA=design_subset$core_RA)

colors_core <- colors_core[colors_core$group %in% df$compartment, ]
df$compartment <- factor(df$compartment, levels=colors_core$group)

shapes_core <- shapes_core[shapes_core$group %in% df$compartment, ]

p <- ggplot(df, aes(x=compartment, y=core_RA, color=compartment, shape=compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(2), alpha=shannon_alpha) +
            scale_y_continuous(labels=scales::percent) +
            scale_color_manual(values=as.character(colors_core$color)) +
            scale_shape_manual(values=shapes_core$shape) +
            labs(x="Compartment", y="Aggregated Relative Abundance") +
            theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
            main_theme

ggsave(paste(figures.dir, "nonflooded_NC_core.pdf", sep=""), p, width=shannon_width, height=shannon_height)

