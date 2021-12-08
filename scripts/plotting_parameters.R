
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

pcoa_width <- 7
pcoa_height <- 6
pcoa_size <- 3
pcoa_alpha <- 0.9
segment_alpha <- 0.4

shannon_width <- 4
shannon_height <- pcoa_height
shannon_alpha <- 0.7

boxplot_size <- 1
boxplot_width <- 0.75
boxplot_jitter_size <- 1.5

width_rec_barplot <- 5
height_rec_barplot <- 3
size_rec_barplot <- 0.35

size_cumsum <- 0.75

proteo_color <- "#2e7f59" 
actino_color <- "#ffff00"
bacteroidetes_color <- "#0000ff"
firmicutes_color <- "#ff4500"
chlorophyta_color <- "#32a852"

at_color <- "#f8766c"
lj_color <- "#00bfc4"
lj_mutant_color <- "#8aaeb6"
cr_color  <- "#32a852"
soil_color <- "#654321"
rhizosphere_color <- c_orange
bacteria_color  <- "#4a80ff"
bacteria_ap_color  <- "#ff8200"
input_color  <- c_black

cr_only_color <- "#16a11d"
cr_only_color <- "#008080"
teal <- "#008080"
grass <- "#536d06"

chlamy_1030_color <- "#3f5e34" 
chlamy_0911_color <- "#5422c9" 
chlamy_0968_color <- "#1c7980"
    
spirogloea_color <- "#eb8f34"
microthamnion_color <- "#c9c547"
klebsormidium_1121_color <- "#3d54a6"

input_shape <- 16
root_shape <- 19
rhizosphere_shape <- 3
soil_shape <- 18
soil_surface_shape <- 25
phycosphere_shape <- 17
aquatic_shape <- 16

covered_color <- "black"
uncovered_color <- "#fcba03"

# ggplot2 theme

main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line=element_line(color="black", size=1),
              axis.ticks=element_line(color="black", size=1),
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(size=18, color="black"),
              legend.position="none")

