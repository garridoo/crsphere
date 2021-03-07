
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# SynCom flask dataset

source("load_data_SC_comp_flowpot.R")

# subset samples of interest

idx <- design$system %in% c("FP", "FPL") &
       # design$condition %in% c("At+Full", "C+Full", "Full", "Full_MS", "Full_TP") &
       design$condition %in% c("At+Full", "At+IRL", "At+ICL", "C+Full", "C+IRL", "C+ICL", "Full_MS", "Full_TP", "Full", "IRL", "ICL") &
       # design$treatment %in% c("Full") &
       design$fraction %in% c("root", "input", "cells", "soil") &
       T

design$fraction_treatment <- paste(design$fraction, design$treatment)

design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]
asv_table_subset_counts <- asv_table_subset_counts[rowSums(asv_table_subset_counts)!=0, ]

# aggregate tables to the family level to compare across SynCom treatments

family_table_norm <- aggregate(asv_table_subset, by=list(taxonomy$family), FUN=sum)
rownames(family_table_norm) <- family_table_norm[, 1]
family_table_norm <- family_table_norm[, -1]

# asv_table_subset <- family_table_norm

### beta diversity

colors <- data.frame(group=c("input", "Cr", "At", "soil"),
                     color=c(input_color, cr_color, at_color, soil_color))

shapes <- data.frame(group=c("root Full", "root ICL", "root IRL", "soil Full", "soil ICL", "soil IRL", "cells Full", "input Full", "input IRL", "input ICL"),
                     shape=c(19, 10, 13, 15, 12, 7, 17, 19, 10, 13))

# PCoA Bray-Curtis

bray_curtis <- vegdist(t(asv_table_subset), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_subset[match(rownames(points), design_subset$SampleID), ])

colors <- colors[colors$group %in% points$host, ]
points$host <- factor(points$host, levels=colors$group)

shapes <- shapes[shapes$group %in% points$fraction_treatment, ]
points$fraction_treatment <- factor(points$fraction_treatment, levels=shapes$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=host, shape=fraction_treatment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle("PCoA of Bray-Curtis distances") +
     main_theme +
     theme(legend.position="right")

ggsave(paste(figures.dir, "competition_flowpot_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### CPCoA analysis

sqrt_transform <- T

capscale.gen <- capscale(bray_curtis ~ host + Condition(treatment * bio_replicate * experiment * system), data=design_subset, add=F, sqrt.dist=sqrt_transform)

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
points_cpcoa <- points_cpcoa[points_cpcoa$compartment!="input", ]

colors_cpcoa <- colors[colors$group %in% points_cpcoa$host, ]
points_cpcoa$host <- factor(points_cpcoa$host, levels=colors$group)

shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$fraction_treatment, ]
points_cpcoa$fraction_treatment <- factor(points_cpcoa$fraction_treatment, levels=shapes$group)

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ host, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('host','seg_x','seg_y')), by='host', sort=FALSE)


# plot CPCo 1 and 2

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=host, shape=fraction_treatment)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "competition_flowpot_CPCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### shannon index

asv_table_subset_rarefied <- rrarefy(asv_table_subset_counts, sample=500)
index <- diversity(t(asv_table_subset_rarefied), index="shannon")
 
index <- cbind(index, design_subset[match(names(index), design_subset$SampleID), ])

colors_index <- colors[colors$group %in% index$host, ]
index$host <- factor(index$host, levels=colors$group)

shapes_index <- shapes[shapes$group %in% index$fraction, ]
index$fraction <- factor(index$fraction, levels=shapes$group)

p <- ggplot(index, aes(x=compartment, y=index, color=host, shape=fraction)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(2), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0, max(index$index))) +
            scale_colour_manual(values=as.character(colors_index$color)) +
            scale_shape_manual(values=shapes_index$shape) +
            labs(x="", y="Shannon index") +
            ggtitle("Shannon diversity") +
            main_theme

ggsave(paste(figures.dir, "competition_flowpot_shannon.pdf", sep=""), p, width=shannon_width, height=shannon_height)

