
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# SynCom flask dataset

source("load_data_SC_flask.R")

# subset samples of interest

idx <- T

design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]
asv_table_subset_counts <- asv_table_subset_counts[rowSums(asv_table_subset_counts)!=0, ]

### beta diversity

colors <- data.frame(group=c("SC", "SC+C", "Input"),
                     color=c(bacteria_color, cr_color, input_color))

shapes <- data.frame(group=c("Day0", "Day1", "Day4", "Day7"),
                     shape=c(18, 3, 19, 17))

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

shapes <- shapes[shapes$group %in% points$timepoint, ]
points$timepoint <- factor(points$timepoint, levels=shapes$group)

# plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=condition, shape=timepoint)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle("PCoA of Bray-Curtis distances") +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "flask_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### CPCoA analysis

sqrt_transform <- T

capscale.gen <- capscale(bray_curtis ~ condition + Condition(bio_replicate * tech_replicate), data=design_subset, add=F, sqrt.dist=sqrt_transform)

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

colors_cpcoa <- colors[colors$group %in% points_cpcoa$condition, ]
points_cpcoa$condition <- factor(points_cpcoa$condition, levels=colors$group)

shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$timepoint, ]
points_cpcoa$timepoint <- factor(points_cpcoa$timepoint, levels=shapes$group)

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ condition, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('condition','seg_x','seg_y')), by='condition', sort=FALSE)

# plot CPCo 1 and 2

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=condition, shape=timepoint)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "flask_CPCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### stacked barplots

df <- melt(asv_table_subset)
colnames(df) <- c("StrainID", "SampleID", "RA")

idx <- order(design$condition, design$timepoint, design$bio_replicate, design$tech_replicate)
df$SampleID <- factor(df$SampleID, levels=design$SampleID[idx])

df$taxonomy <- taxonomy$order[match(df$StrainID, taxonomy$ID)]

# plot barplots

p <- ggplot(df, aes(x=SampleID, y=RA, fill=taxonomy)) +
     geom_bar(stat="identity", color="black") +
     scale_y_continuous(limits=c(0, 1), labels=percent) +
     labs(y="Relative abundance") +
     main_theme +
     theme(axis.text.x=element_text(size=10, angle=90, hjust=1)) +
     theme(legend.position="right")

ggsave(paste(figures.dir, "flask_barplot.pdf", sep=""), p, width=16, height=6)

### shannon index

asv_table_subset_rarefied <- rrarefy(asv_table_subset_counts, sample=500)
index <- diversity(t(asv_table_subset_rarefied), index="shannon")
 
index <- cbind(index, design_subset[match(names(index), design_subset$SampleID), ])

colors_index <- colors[colors$group %in% index$condition, ]
index$condition <- factor(index$condition, levels=colors$group)

shapes_index <- shapes[shapes$group %in% index$timepoint, ]
index$timepoint <- factor(index$timepoint, levels=shapes$group)

p <- ggplot(index, aes(x=condition, y=index, color=condition, shape=timepoint)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            geom_jitter(position=position_jitterdodge(2), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0, max(index$index))) +
            scale_colour_manual(values=as.character(colors_index$color)) +
            scale_shape_manual(values=shapes_index$shape) +
            labs(x="", y="Shannon index") +
            ggtitle("Shannon diversity") +
            main_theme

ggsave(paste(figures.dir, "flask_shannon.pdf", sep=""), p, width=shannon_width, height=shannon_height)

