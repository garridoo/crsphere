
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load At Cr natural community dataset

source("load_data_SC_split.R")

# subset samples of interest

idx <- design$Content %in% c("C+SC", "SC", "Input") &
       # design$Content %in% c("C+SC", "Input") &
       # design$Content %in% c("SC", "Input") &
       T

design_subset <- design[idx, ]
asv_table_subset <- asv_table_norm[, idx]
asv_table_subset_counts <- asv_table[, idx]

asv_table_subset_counts <- asv_table_subset_counts[rowSums(asv_table_subset_counts)!=0, ]

mod <- model.matrix(~condition, data=design_subset)
combat_asv_table <- ComBat(dat=as.matrix(asv_table_subset_counts), batch=design_subset$split_exp, mod=mod, par.prior=F, mean.only=T, ref.batch=NULL)
idx <- rowSums(is.na(combat_asv_table))==0
combat_asv_table <- combat_asv_table[idx, ]

asv_table_counts <- round(combat_asv_table * 1000)
asv_table_subset <- apply(asv_table_counts, 2, function(x) x/sum(x))

### beta diversity

# define colors (by content) and shapes (by neighbor)

colors <- data.frame(group=c("SC/C", "C+SC/TP", "C+SC/SC", "SC/C+SC","C+SC/C", "SC/TP", "Input/"),
                     color=c(bacteria_color, cr_color, cr_color, bacteria_color, cr_color, bacteria_color, input_color))

shapes <- data.frame(group=c("SC/C", "C+SC/TP", "C+SC/SC", "SC/C+SC","C+SC/C", "SC/TP", "Input/"),
                     shape=c(16, 1, 17, 15, 16, 1, 16))

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
points$community <- factor(points$condition, levels=colors$group)

shapes <- shapes[shapes$group %in% points$condition, ]
points$condition <- factor(points$condition, levels=shapes$group)
 
# # plot PCo 1 and 2

p <- ggplot(points, aes(x=x, y=y, color=condition, shape=condition)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     scale_colour_manual(values=as.character(colors$color)) +
     scale_shape_manual(values=shapes$shape) +
     labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle("PCoA of Bray-Curtis distances") +
     main_theme +
     theme(legend.position="right")

ggsave(paste(figures.dir, "split_PCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### CPCoA analysis

sqrt_transform <- T

if ("SC" %in% design_subset$Content) {

    # formula for all conditions or SC content
    capscale.gen <- capscale(bray_curtis ~ chlamy + Condition(split_exp * Replicate), data=design_subset, add=F, sqrt.dist=sqrt_transform)

}

if (!"SC" %in% design_subset$Content) {

    # formula for SC+C content
    capscale.gen <- capscale(bray_curtis ~ condition + Condition(split_exp * Replicate), data=design_subset, add=F, sqrt.dist=sqrt_transform)

}

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

shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$condition, ]
points_cpcoa$condition <- factor(points_cpcoa$condition, levels=shapes$group)

# calculate centroids per condition and joining segments

centroids <- aggregate(cbind(points_cpcoa$x, points_cpcoa$y) ~ condition, data=points_cpcoa, FUN=mean)
segments <- merge(points_cpcoa, setNames(centroids, c('condition','seg_x','seg_y')), by='condition', sort=FALSE)

# plot CPCo 1 and 2

p <- ggplot(points_cpcoa, aes(x=x, y=y, color=condition, shape=condition, label=SampleID)) +
     geom_point(alpha=pcoa_alpha, size=pcoa_size) +
     geom_segment(data=segments, mapping=aes(xend=seg_x, yend=seg_y), alpha=segment_alpha) +
     # geom_text() + 
     scale_colour_manual(values=as.character(colors_cpcoa$color)) +
     scale_shape_manual(values=shapes_cpcoa$shape) +
     labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
          y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
     ggtitle(paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
     main_theme +
     theme(legend.position="none")

ggsave(paste(figures.dir, "split_CPCoA.pdf", sep=""), p, width=pcoa_width, height=pcoa_height)

### shannon index

asv_table_subset_rarefied <- rrarefy(asv_table_subset_counts, sample=500)
index <- diversity(t(asv_table_subset_rarefied), index="shannon")
 
index <- cbind(index, design_subset[match(names(index), design_subset$SampleID), ])

colors_index <- colors[colors$group %in% index$condition, ]
index$condition <- factor(index$condition, levels=colors$group)

shapes_index <- shapes[shapes$group %in% index$condition, ]
index$condition <- factor(index$condition, levels=shapes$group)

p <- ggplot(index, aes(x=condition, y=index, color=condition, shape=condition)) +
            geom_boxplot(alpha=1, outlier.size=0, size=boxplot_size, width=boxplot_width, fill="transparent") +
            stat_summary(fun.y="mean",geom="crossbar") +
            geom_jitter(position=position_jitterdodge(2), size=boxplot_jitter_size, alpha=shannon_alpha) +
            scale_y_continuous(limits=c(0, max(index$index))) +
            scale_colour_manual(values=as.character(colors_index$color)) +
            scale_shape_manual(values=shapes_index$shape) +
            labs(x="", y="Shannon index") +
            ggtitle("Shannon diversity") +
            main_theme

ggsave(paste(figures.dir, "split_shannon.pdf", sep=""), p, width=shannon_width, height=shannon_height)

