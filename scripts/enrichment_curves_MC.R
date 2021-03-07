
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# load data and calculate enriched ASVs

source("enrichment_MC.R")

# ancillary functions

sem <- function(x) {sd(x)/length(x)}

# subset samples

idx <- design$system=="flask" &
       design$community %in% c("Input", "B+C", "B") &
       design$soil %in% c("CAS") &
       design$experiment %in% c("AL1", "AL2", "AL3", "AL4", "AL5", "AL6") &
       design$medium %in% c("TP", "TP+AP", "Input") &
       TRUE

design_subset <- design[idx, ]
asv_table_subset <- asv_table[, idx]
asv_table_subset_norm <- asv_table_norm[, idx]

# define colors, shapes, and sizes

colors <- data.frame(group=c("RA_Input_ASVs", "RA_BC_ASVs", "RA_B_ASVs", "RA_BAP_ASVs"),
                     color=c(soil_color, cr_color, bacteria_color, bacteria_ap_color))

sizes <- data.frame(group=c("Day0", "Day1", "Day4", "Day7", "Day8", "Day11"),
                    size=c(1.50, 1.75, 2.00, 2.25, 2.50, 2.75))

# remove empty columns from the ASV table

idx <- rowSums(asv_table_subset_norm)!=0
asv_table_subset <- asv_table_subset[idx, ]
asv_table_subset_norm <- asv_table_subset_norm[idx, ]

# create data frame of RAs

df <- melt(asv_table_subset_norm)
colnames(df) <- c("ASV", "SampleID", "RA")
df <- cbind(df, design[match(df$SampleID, design$SampleID), -1])

# add information about each ASV enrichment status and aggregate RAs

df$enriched_Input <- FALSE
df$enriched_Input[df$ASV %in% Input_ASVs] <- TRUE
df_aggregated <- aggregate(df$RA, by=list(df$SampleID, df$enriched_Input), FUN=sum)
colnames(df_aggregated) <- c("SampleID", "enrichment", "RA_Input_ASVs")

df$enriched_BC <- FALSE
df$enriched_BC[df$ASV %in% BC_ASVs] <- TRUE
df_aggregated$RA_BC_ASVs <- aggregate(df$RA, by=list(df$SampleID, df$enriched_BC), FUN=sum)[, 3]

df$enriched_B <- FALSE
df$enriched_B[df$ASV %in% B_ASVs] <- TRUE
df_aggregated$RA_B_ASVs <- aggregate(df$RA, by=list(df$SampleID, df$enriched_B), FUN=sum)[, 3]

df$enriched_BAP <- FALSE
df$enriched_BAP[df$ASV %in% BAP_ASVs] <- TRUE
df_aggregated$RA_BAP_ASVs <- aggregate(df$RA, by=list(df$SampleID, df$enriched_BAP), FUN=sum)[, 3]

df_aggregated <- df_aggregated[df_aggregated$enrichment, -2]
df_aggregated <- melt(df_aggregated)
colnames(df_aggregated) <- c("SampleID", "enrichment", "RA")

# create data frame of aggregated RAs (by enrichment)

df_aggregated <- cbind(df_aggregated, design_subset[match(df_aggregated$SampleID, design_subset$SampleID), -1])

# subset B+C samples

df_cond <- df_aggregated[df_aggregated$community %in% c("Input", "B+C") & df_aggregated$medium %in% c("Input", "TP", "BD"), ]

# aggregate and calculate statistics

df_stats_cond <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=mean)
colnames(df_stats_cond) <- c("enrichment", "timepoint", "mean")
df_stats_cond$sd <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=sd)[, 3]

# adjust colors

colors_BC  <- colors[colors$group %in% df_stats_cond$enrichment, ]
df_stats_cond$colors_BC <- factor(df_stats_cond$enrichment, level=colors_BC$group)

sizes_BC <- sizes[sizes$group %in% df_stats_cond$timepoint, ]
df_stats_cond$sizes_BC <- factor(df_stats_cond$timepoint, level=sizes_BC$group)

df_stats_cond$timepoint <- factor(df_stats_cond$timepoint, levels=sizes_BC$group)

# plot curves

p1 <- ggplot(df_stats_cond, aes(x=timepoint, y=mean, color=enrichment)) +
            geom_line(aes(group=enrichment), alpha=0.7, size=1.5) +
            geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, fatten=1.5), alpha=0.7, size=1.5) +
            scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Condition", y="") +
            ggtitle("B+C") +
            main_theme +
            theme(legend.position="right")

# subset B samples

df_cond <- df_aggregated[df_aggregated$community %in% c("Input", "B") & df_aggregated$medium %in% c("Input", "TP", "BD"), ]

# aggregate and calculate statistics

df_stats_cond <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=mean)
colnames(df_stats_cond) <- c("enrichment", "timepoint", "mean")
df_stats_cond$sd <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=sd)[, 3]

# adjust colors

colors_B  <- colors[colors$group %in% df_stats_cond$enrichment, ]
df_stats_cond$colors_B <- factor(df_stats_cond$enrichment, level=colors_B$group)

sizes_B <- sizes[sizes$group %in% df_stats_cond$timepoint, ]
df_stats_cond$sizes_B <- factor(df_stats_cond$timepoint, level=sizes_B$group)

df_stats_cond$timepoint <- factor(df_stats_cond$timepoint, levels=sizes_B$group)

# plot curves

p2 <- ggplot(df_stats_cond, aes(x=timepoint, y=mean, color=enrichment)) +
            geom_line(aes(group=enrichment), alpha=0.7, size=1.5) +
            geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, fatten=1.5), alpha=0.7, size=1.5) +
            scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Condition", y="") +
            ggtitle("B") +
            main_theme +
            theme(legend.position="right")

# subset B (+AP) samples

df_cond <- df_aggregated[df_aggregated$community %in% c("Input", "B") & df_aggregated$medium %in% c("Input", "TP+AP", "BD+AP"), ]

# aggregate and calculate statistics

df_stats_cond <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=mean)
colnames(df_stats_cond) <- c("enrichment", "timepoint", "mean")
df_stats_cond$sd <- aggregate(df_cond$RA, by=list(df_cond$enrichment, df_cond$timepoint), FUN=sd)[, 3]

# adjust colors

colors_BAP  <- colors[colors$group %in% df_stats_cond$enrichment, ]
df_stats_cond$colors_BAP <- factor(df_stats_cond$enrichment, level=colors_BAP$group)

sizes_BAP <- sizes[sizes$group %in% df_stats_cond$timepoint, ]
df_stats_cond$sizes_BAP <- factor(df_stats_cond$timepoint, level=sizes_BAP$group)

df_stats_cond$timepoint <- factor(df_stats_cond$timepoint, levels=sizes_BAP$group)

# plot curves

p3 <- ggplot(df_stats_cond, aes(x=timepoint, y=mean, color=enrichment)) +
            geom_line(aes(group=enrichment), alpha=0.7, size=1.5) +
            geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd, fatten=1.5), alpha=0.7, size=1.5) +
            scale_y_continuous(limits=c(0, 1), labels=scales::percent) +
            scale_colour_manual(values=as.character(colors$color)) +
            labs(x="Condition", y="") +
            ggtitle("BAP") +
            main_theme +
            theme(legend.position="right")

pg1 <- arrangeGrob(p1, p2, p3, ncol=1, nrow=3, heights=c(4, 4, 4), widths=12)

ggsave(paste(figures.dir, "enrichment_curves_MC.pdf", sep=""), pg1, width=12, height=12)

