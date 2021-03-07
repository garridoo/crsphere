
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load At Cr natural community dataset

source("load_data_NC.R")

# subset samples of interest

idx <- design$system=="soil" &
       (design$compartment %in% c("root", "rhizosphere", "phycosphere", "soil") & design$timepoint=="Day36") &
       T
design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]
taxonomy_subset <- taxonomy[idx, ]

### beta diversity

colors <- data.frame(group=c("B", "B+At", "B+C"),
                     color=c(soil_color, at_color, cr_color))

shapes <- data.frame(group=c("root", "rhizosphere", "phycosphere", "soil"),
                     shape=c(root_shape, rhizosphere_shape, phycosphere_shape, soil_shape))

# PCoA Bray-Curtis

bray_curtis <- vegdist(t(asv_table_subset), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_subset[match(rownames(points), design_subset$SampleID), ])

colors <- colors[colors$group %in% points$community, ]
points$community <- factor(points$community, levels=colors$group)

shapes <- shapes[shapes$group %in% points$compartment, ]
points$compartment <- factor(points$compartment, levels=shapes$group)

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points$x, points$y) ~ compartment, data=points, FUN=mean)
segments <- merge(points, setNames(centroids, c('compartment','seg_x','seg_y')), by='compartment', sort=FALSE)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=community, shape=compartment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle("PCoA of Bray-Curtis distances") +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "AtCr_NC_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### CPCoA analysis

design_subset$host_associated <- TRUE
design_subset$host_associated[design_subset$compartment %in% c("soil", "rhizosphere")] <- FALSE

sqrt_transform <- T

capscale.gen <- capscale(bray_curtis ~ compartment + Condition(replicate * timepoint),
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

colors_cpcoa <- colors[colors$group %in% points_cpcoa$community, ]
points_cpcoa$community <- factor(points_cpcoa$community, levels=colors$group)

shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$compartment, ]
points_cpcoa$compartment <- factor(points_cpcoa$compartment, levels=shapes$group)

# plot CPCo 1 and 2

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ compartment, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('compartment','seg_x','seg_y')), by='compartment', sort=FALSE)

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=community, shape=compartment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "AtCr_NC_CPCoA_1_2.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

# plot CPCo 2 and 3

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$y, points_cpcoa$z) ~ compartment, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('compartment','seg_x','seg_y')), by='compartment', sort=FALSE)

p <- ggplot(points_cpcoa, aes(x=y, y=z, color=community, shape=compartment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "AtCr_NC_CPCoA_2_3.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### shannon index

asv_table_subset_rarefied <- rrarefy(asv_table_subset_counts, sample=500)
index <- diversity(t(asv_table_subset_rarefied), index="shannon")
 
index <- cbind(index, design_subset[match(names(index), design_subset$SampleID), ])

colors_index <- colors[colors$group %in% index$community, ]
index$community <- factor(index$community, levels=colors$group)

shapes_index <- shapes[shapes$group %in% index$compartment, ]
index$compartment <- factor(index$compartment, levels=shapes$group)

# reorder boxplots

index$compartment <- factor(index$compartment, levels=rev(shapes$group))
shapes_index$shape  <- rev(shapes_index$shape)

p <- ggplot(index, aes(x=compartment, y=index, color=community, shape=compartment)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(2), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0, max(index$index))) +
            scale_colour_manual(values=as.character(colors_index$color)) +
            scale_shape_manual(values=shapes_index$shape) +
            labs(x="Compartment", y="Shannon index") +
            ggtitle("Shannon diversity") +
            main_theme

ggsave(paste(figures.dir, "AtCr_NC_shannon.pdf", sep=""), p, width=shannon_width, height=shannon_height)

