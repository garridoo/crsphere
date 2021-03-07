 
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

source("load_data_MC.R")

# subset samples of interest

idx <- design$system=="flask" &
       design$community %in% c("B+C", "B", "B+At", "Input") &
       design$experiment %in% c("AL1", "AL2", "AL3", "AL4", "AL5", "AL6") &
       ((design$medium %in% c("TP", "BD") & design$community == "B+C") | (design$medium %in% c("TP", "BD", "TP+AP", "BD+AP") & design$community == "B") | design$medium == "Input" | design$medium == "agar" ) &
       TRUE

design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]
taxonomy_subset <- taxonomy[idx, ]

design_subset$condition  <- paste(design_subset$community, design_subset$AP)

### beta diversity

colors <- data.frame(group=c("B+C FALSE", "B FALSE", "B TRUE", "Input FALSE"),
                     color=c(cr_color, bacteria_color, bacteria_ap_color, soil_color))

shapes <- data.frame(group=c("CAS", "Golm"),
                     shape=c(19, 17))

sizes <- data.frame(group=c("Day0", "Day1", "Day4", "Day7", "Day8", "Day11"),
                    size=c(1.50, 1.75, 2.00, 2.25, 2.50, 2.75))

# PCoA Bray-Curtis

bray_curtis <- vegdist(t(asv_table_subset), method="bray")

k <- 2
pcoa <- cmdscale(bray_curtis, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig
points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points <- cbind(points, design_subset[match(rownames(points), design_subset$SampleID), ])

colors <- colors[colors$group %in% points$condition, ]
points$condition <- factor(points$condition, levels=colors$group)

shapes <- shapes[shapes$group %in% points$soil,  ]
points$soil <- factor(points$soil, levels=shapes$group)

sizes  <- sizes[sizes$group %in% points$timepoint, ]
points$timepoint  <- factor(points$timepoint, level=sizes$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=condition, shape=soil, size=timepoint)) +
     geom_point(alpha=pcoa_alpha) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     scale_size_manual(values=sizes$size) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle("PCoA of Bray-Curtis distances") +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "AtCr_MC_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### CPCoA analysis

sqrt_transform <- T

capscale.gen <- capscale(bray_curtis ~ community * medium + Condition(experiment * timepoint * soil), data=design_subset, add=F, sqrt.dist=sqrt_transform)

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

colors <- colors[colors$group %in% points_cpcoa$condition, ]
points_cpcoa$condition <- factor(points_cpcoa$condition, levels=colors$group)

shapes <- shapes[shapes$group %in% points_cpcoa$soil,  ]
points_cpcoa$soil <- factor(points_cpcoa$soil, levels=shapes$group)

sizes  <- sizes[sizes$group %in% points_cpcoa$timepoint, ]
points_cpcoa$timepoint  <- factor(points_cpcoa$timepoint, level=sizes$group)

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ condition, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('condition','seg_x','seg_y')), by='condition', sort=FALSE)

# plot CPCo 1 and 2

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=condition, shape=soil, size=timepoint)) +
     geom_point(alpha=pcoa_alpha) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha, size=.5) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     scale_size_manual(values=sizes$size) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "AtCr_MC_CPCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### shannon index

asv_table_subset_rarefied <- rrarefy(asv_table_subset_counts, sample=500)
index <- diversity(t(asv_table_subset_rarefied), index="shannon")
 
index <- cbind(index, design_subset[match(names(index), design_subset$SampleID), ])

colors <- colors[colors$group %in% index$condition, ]
index$condition <- factor(index$condition, levels=colors$group)

shapes <- shapes[shapes$group %in% index$soil,  ]
index$soil <- factor(index$soil, levels=shapes$group)

sizes <- sizes[sizes$group %in% index$timepoint, ]
index$sizes <- factor(index$timepoint, level=sizes$group)

index$timepoint <- factor(index$timepoint, levels=sizes$group)

# reorder boxplots

p <- ggplot(index, aes(x=timepoint, y=index, color=condition, shape=soil)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(0.35), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0, max(index$index))) +
            scale_colour_manual(values=as.character(colors$color)) +
            scale_shape_manual(values=shapes$shape) +
            labs(x="Condition", y="Shannon index") +
            ggtitle("Shannon diversity") +
            main_theme +
            theme(legend.position="right")

ggsave(paste(figures.dir, "AtCr_MC_shannon.pdf", sep=""), p, width=12, height=shannon_height)

