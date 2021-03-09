## Scripts of Duran, Flores-Uribe, *et al.*, Characterization of the *Chlamydomonas reinhardtii* phycosphere reveals conserved features of the plant microbiota

Correspondence to [Ruben Garrido-Oter](mailto:garridoo@mpipz.mpg.de).

These scripts are made available to facilitate the reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper. Raw data and intermediate results necessary to run these scripts can also be downloaded [here](http://www.mpipz.mpg.de/scripts).

Our manuscript is currently available as a preprint in [*bioRxiv*](https://doi.org/10.1101/2021.03.04.433956).

---------------------------

To run the scripts contained in this repository, we recommend to download first the raw and intermediate data from [here](http://www.at-sphere.com/cr.tar.gz). After extracting its contents, the script [paths.R](https://github.com/garridoo/crsphere/blob/master/scripts/paths.R) can be modified so that the base path points to the directory containing the data. After installing the necessary R libraries, all experiments can be analyzed independently by loading the corresponding data and subsequently running the remastering scripts.

### Accession numbers

Raw *16S* rRNA amplicon and whole-genome sequencing reads were deposited in the European Nucleotide Archive (ENA) under the accession number [PRJEB43117](XXX).

In addition, raw sequencing data from all SynCom experiments (FASTQ files and mapping files including barcodes, biological and technical metadata, etc.) along with intermediate results (ASV tables, etc.) can be downloaded [here](http://www.at-sphere.com/cr.tar.gz).

### Analyses of culture-independent *16S* rRNA amplicon data

#### Greenhouse experiment

The scripts [load_data_NC.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_NC.R),
[diversity_NC.R](https://github.com/garridoo/crsphere/blob/master/scripts/diversity_NC.R),
[enrichment_NC.R](https://github.com/garridoo/crsphere/blob/master/scripts/enrichment_NC.R), and
[enrichment_curves_NC.R](https://github.com/garridoo/crsphere/blob/master/scripts/enrichment_curves_NC.R) are used to load data from the natural community (greenhouse experiment; Fig. S1A) and perform alpha- and beta-diversity analyses and identify AVSs enriched in the different compartments. These scripts correspond to the analyses shown in Fig. 1, and Fig. S2 of the manuscript.

#### Root and phycosphere meta-analysis

The scripts [load_data_yeoh.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_yeoh.R),
[load_data_thiergart_harbort.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_thiergart_harbort.R),
[core.R](https://github.com/garridoo/crsphere/blob/master/scripts/core.R), and
[core_plots.R](https://github.com/garridoo/crsphere/blob/master/scripts/core_plots.R) are used to load the *16S* rRNA amplicon datasets published in Yeoh *et al.*, 2017, Thiergart *et al.*, 2019, and Harbort *et al.*, 2020, and perform the analyses of core microbial taxa described in the manuscript and depicted in Fig. 2.

#### Mesocosm experiments

To analyse the data from the mesocosm experiments (Fig. S1B), including alpha- and beta-diversity, enrichment tests and assessment of algal growth data, the following scripts can be used: [load_data_MC.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_MC.R),
[diversity_MC.R](https://github.com/garridoo/crsphere/blob/master/scripts/diversity_MC.R),
[enrichment_MC.R](https://github.com/garridoo/crsphere/blob/master/scripts/enrichment_MC,R),
[enrichment_curves_MC.R](https://github.com/garridoo/crsphere/blob/master/scripts/enrichment_curves_MC.R), and
[cell_counts_MC.R](https://github.com/garridoo/crsphere/blob/master/scripts/cell_counts_MC.R). These scripts generate the analyses shown in Fig. 3.

#### Sequence-Indexed Phycobacterial Library (IPL)

The scripts [IPL_load_data.R](https://github.com/garridoo/crsphere/blob/master/scripts/IPL_load_data.R),
[IPL.R](https://github.com/garridoo/crsphere/blob/master/scripts/IPL.R), and
[IPL_recovery_rates.R](https://github.com/garridoo/crsphere/blob/master/scripts/IPL_recovery_rates.R) are used to used to cross-reference data from shallow sequencing of the IPL (Fig S1C) and culture-independent profiling of start inocula (e.g., greenhouse data), and to estimate species recovery rates (Fig. 1D, and Fig. S4A).

### Analyses of culture-dependent *16S* rRNA amplicon data

#### SynCom flask experiment

The scripts [load_data_SC_flask.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_SC_flask.R),
[diversity_SC_flask.R](https://github.com/garridoo/crsphere/blob/master/scripts/diversity_SC_flask.R), and
[fluorescence_SC_flask.R](https://github.com/garridoo/crsphere/blob/master/scripts/fluorescence_SC_flask.R) can be used to load data from the flask SynCom experiment (Fig. SD), and perform the diversity analyses shown in Fig. 4 of the manuscript.

#### SynCom FlowPot and flask competition experiments

To process the data from the competition experiments in the soil-based gnotobiotic system (Fig. S1E), the following scripts can be used:  [load_data_SC_comp_flowpot.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_SC_comp_flowpot.R),
[diversity_SC_comp_flowpot.R](https://github.com/garridoo/crsphere/blob/master/scripts/diversity_SC_comp_flowpot.R), and
[boxplots_SC_comp_flowpot.R](https://github.com/garridoo/crsphere/blob/master/scripts/boxplots_SC_comp_flowpot.R).
Likewise, the data from the competition experiment in a liquid based system can be processes with the scripts [load_data_SC_comp_flask.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_SC_flask.R),
[boxplots_SC_comp_flask.R](https://github.com/garridoo/crsphere/blob/master/scripts/boxplots_SC_comp_flask.R), and
[cell_counts_SC_flask.R](https://github.com/garridoo/crsphere/blob/master/scripts/cell_counts_SC_flask.R). Together, these scripts generate the panels shown in Fig. 5.

#### SynCom split system experiment

The scripts [load_data_SC_split.R](https://github.com/garridoo/crsphere/blob/master/scripts/load_data_SC_split.R),
[diversity_SC_split.R](https://github.com/garridoo/crsphere/blob/master/scripts/diversity_SC_split.R),
[fluorescence_SC_split.R](https://github.com/garridoo/crsphere/blob/master/scripts/fluorescence_SC_split.R), and
[cell_counts_SC_split.R](https://github.com/garridoo/crsphere/blob/master/scripts/cell_counts_SC_split.R) can be used to process the data from the SynCom split system (Fig. S1F) and perform the analyses depicted in Fig. 6 of the manuscript.

### Whole-genome assembly, quality control, and annotation of the *Cr*-SPHERE core culture collection

Raw sequencing data (FASTQ files) from the core *Cr*-SPHERE culture collection as well as assemblies (FNA files), nucleotide ORFs (FFN), amino acid ORFs (FAA), GFF files, KEGG annotations (KO), reference *16S* rRNA sequences, AMPHORA marker gene alignments, and metadata can be downloaded in bulk [here](http://www.at-sphere.com/cr.tar.gz). The [scripts](https://github.com/garridoo/ljsphere) used in Wippel *et al.*, 2020 were repurposed for these analyses and are described below.

[assembly.functions.sh](https://github.com/garridoo/ljsphere/blob/master/assembly.functions.sh): bash script containing auxiliary functions for whole-genome assembly.

[assembly.sh](https://github.com/garridoo/ljsphere/blob/master/assembly.sh): script used to assemble the genomes using SOAPdenovo and A5. It can be run in parallel using either the custom script bellow or the gnu parallel suit.

[parallel.sh](https://github.com/garridoo/lsphere/blob/master/parallel.sh): custom script to run bash functions in parallel in a multi-core machine.

[assembly_stats.R](https://github.com/garridoo/ljsphere/blob/master/assembly_stats.R): R script used to generate assembly statistics as well as GC and *k*-mer spectral projections. The output of this script contains clean assemblies (all contigs < 1,000 bp are removed) as well as a PDF file containing a report which was used to manually inspect for likely contaminated assemblies.

### Additional scripts

[plotting_parameters.R](https://github.com/garridoo/crsphere/blob/master/scripts/plotting_parameters.R), and [plotting_functions.R](https://github.com/garridoo/crsphere/blob/master/scripts/plotting_functions.R) contain R scripts containing plotting parameters such as colors, ggplot2 themes, etc. as well as auxiliary functions for generating the figures reported in the paper (some of the details may vary with respect to the published version).

---------------------------

For any questions regarding these scripts, please contact

Ruben Garrido-Oter

garridoo@mpipz.mpg.de
