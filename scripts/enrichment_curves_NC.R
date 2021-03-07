
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# load data and calculate enriched ASVs

source("enrichment_NC.R")

# ancillary functions

sem <- function(x) {sd(x)/length(x)}

# subset samples

idx <- design$system=="soil" &
       design$compartment %in% c("soil", "phycosphere", "root") &
       TRUE

design_subset <- design[idx, ]
asv_table_subset <- asv_table[, idx]
asv_table_subset_norm <- asv_table_norm[, idx]

# define colors, shapes, and sizes

colors <- data.frame(group=c("RA_soil_ASVs", "RA_phycosphere_ASVs", "RA_root_ASVs"),
                     color=c(soil_color, cr_color, at_color))

# remove empty columns from the ASV table

idx <- rowSums(asv_table_subset_norm)!=0
asv_table_subset <- asv_table_subset[idx, ]
asv_table_subset_norm <- asv_table_subset_norm[idx, ]

# create data frame of RAs

df <- melt(asv_table_subset_norm)
colnames(df) <- c("ASV", "SampleID", "RA")
df <- cbind(df, design[match(df$SampleID, design$SampleID), -1])

# add information about each ASV enrichment status and aggregate RAs

df$enriched_soil <- FALSE
df$enriched_soil[df$ASV %in% soil_ASVs] <- TRUE
df_aggregated <- aggregate(df$RA, by=list(df$SampleID, df$enriched_soil), FUN=sum)
colnames(df_aggregated) <- c("SampleID", "enrichment", "RA_soil_ASVs")

df$enriched_phycosphere <- FALSE
df$enriched_phycosphere[df$ASV %in% phycosphere_ASVs] <- TRUE
df_aggregated$RA_phycosphere_ASVs <- aggregate(df$RA, by=list(df$SampleID, df$enriched_phycosphere), FUN=sum)[, 3]

df$enriched_root <- FALSE
df$enriched_root[df$ASV %in% root_ASVs] <- TRUE
df_aggregated$RA_root_ASVs <- aggregate(df$RA, by=list(df$SampleID, df$enriched_root), FUN=sum)[, 3]

df_aggregated <- df_aggregated[df_aggregated$enrichment, -2]
df_aggregated <- melt(df_aggregated)
colnames(df_aggregated) <- c("SampleID", "enrichment", "RA")

# create data frame of aggregated RAs (by enrichment)

df_aggregated <- cbind(df_aggregated, design_subset[match(df_aggregated$SampleID, design_subset$SampleID), -1])

# subset phycosphere samples

df_cond <- df_aggregated[df_aggregated$compartment %in% c("phycosphere") | df_aggregated$timepoint %in% c("Day0"), ]

# aggregate and calculate statistics

df_stats_cond <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=mean)
colnames(df_stats_cond) <- c("enrichment", "timepoint", "mean")
df_stats_cond$sd <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=sd)[, 3]

# adjust colors

colors_phycosphere  <- colors[colors$group %in% df_stats_cond$enrichment, ]
df_stats_cond$colors_phycosphere <- factor(df_stats_cond$enrichment, level=colors_phycosphere$group)

df_stats_cond$timepoint <- factor(df_stats_cond$timepoint, levels=c("Day0", "Day7", "Day14", "Day21", "Day28", "Day36"))

# plot curves

p1 <- ggplot(df_stats_cond, aes(x=timepoint, y=mean, color=enrichment)) +
            geom_line(aes(group=enrichment), alpha=0.7, size=1.5) +
            geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, fatten=1.5), alpha=0.7, size=1.5) +
            scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Condition", y="") +
            ggtitle("phycosphere") +
            main_theme +
            theme(legend.position="right")

# subset root samples

df_cond <- df_aggregated[df_aggregated$compartment %in% c("root") | df_aggregated$timepoint %in% c("Day0"), ]

# aggregate and calculate statistics

df_stats_cond <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=mean)
colnames(df_stats_cond) <- c("enrichment", "timepoint", "mean")
df_stats_cond$sd <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=sd)[, 3]

# adjust colors

colors_root  <- colors[colors$group %in% df_stats_cond$enrichment, ]
df_stats_cond$colors_root <- factor(df_stats_cond$enrichment, level=colors_root$group)

df_stats_cond$timepoint <- factor(df_stats_cond$timepoint, levels=c("Day0", "Day7", "Day14", "Day21", "Day28", "Day36"))

# plot curves

p2 <- ggplot(df_stats_cond, aes(x=timepoint, y=mean, color=enrichment)) +
            geom_line(aes(group=enrichment), alpha=0.7, size=1.5) +
            geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, fatten=1.5), alpha=0.7, size=1.5) +
            scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Condition", y="") +
            ggtitle("root") +
            main_theme +
            theme(legend.position="right")

# subset soil samples

df_cond <- df_aggregated[df_aggregated$compartment %in% c("soil"), ]

# aggregate and calculate statistics

df_stats_cond <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=mean)
colnames(df_stats_cond) <- c("enrichment", "timepoint", "mean")
df_stats_cond$sd <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=sd)[, 3]

# adjust colors

colors_soil  <- colors[colors$group %in% df_stats_cond$enrichment, ]
df_stats_cond$colors_soil <- factor(df_stats_cond$enrichment, level=colors_soil$group)

df_stats_cond$timepoint <- factor(df_stats_cond$timepoint, levels=c("Day0", "Day7", "Day14", "Day21", "Day28", "Day36"))

# plot curves

p3 <- ggplot(df_stats_cond, aes(x=timepoint, y=mean, color=enrichment)) +
            geom_line(aes(group=enrichment), alpha=0.7, size=1.5) +
            geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, fatten=1.5), alpha=0.7, size=1.5) +
            scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Condition", y="") +
            ggtitle("soil") +
            main_theme +
            theme(legend.position="right")

pg1 <- arrangeGrob(p1, p2, p3, ncol=1, nrow=3, heights=c(4, 4, 4), widths=12)

ggsave(paste(figures.dir, "enrichment_curves_NC.pdf", sep=""), pg1, width=12, height=12)

