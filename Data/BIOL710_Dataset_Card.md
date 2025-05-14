**Dataset Title:** Population Variance of Thermoflexus hugenholtzii in
Great Boiling Spring

**Dataset Curator:** Gillyan Valencia, MS graduate student, José R. de
la Torre lab, Department of Biology, San Francisco State University
(SFSU)

**Dataset version:** 1.2, 4.24.25

**Dataset citation and DOI:** 
Data statement author: Gillyan Valencia, MS graduate student, de la Torre Lab, Department of Biology, San Francisco State University

**Data statement citation:** authors, date, title, version, institution, and URL or DOI.

### Languages

English

## Datset Summary

This dataset contains frequencies of single nucleotide variant (SNV) and
single codon variant (SCV) data for selected core genes (rpoB, gyrB,
recA) across multiple environmental samples. Each variant record
includes positional information, associated gene annotations, allele or
codon frequencies, and contextual data such as sampling site, month, and
temperature. This dataset provides a resource for exploring microbial
population structure, genetic diversity, and evolutionary pressures
across spatial and temporal gradients.

## Data Format

All data are frequency values 

## Data Collection:

**Field Sampling:** Data was collected from field sampling of the
terrestrial hot spring Great Boiling Spring in Gerlach, Nevada, United
States. Sediment was collected from the top 1 cm profile of each sample
site into a sterile plastic beaker. The sediment was homogenized and
alloquoted in 5 mL increments into 50 mL sterile falcon tubes. The
samples were stored on dry ice until transportation back to the
laboratory where they were kept at -80 C until DNA extraction. DNA was
extracted from samples and short-read sequencing was obtained through
illumina sequencing.

**Bioinformatics:** Raw reads were obtained of the forward and reverse
strand for each of the 5 sample points through illumina sequencing. The
raw reads were de novo assembled separately using SPAdes (v. 4.1.0) with 
parameters: -meta -k 21,33,55,77,99. The generated scaffolds.fasta file was 
used for binning with MetaBat2 (v. 2.17). The metagenome-assembled 
genomes (MAGs) generated from MetaBat2 were checked for contamination and
completeness to ensure genomes selected were of high quality as well as 
taxonomically classified with gtdbtk (v2.4.0). Thermoflexus was chosen because 
it was ubiquitously dispersed across all sample sites and had high completeness 
and low contamination.

The reference MAG of Thermoflexus hugenholtzii was chosen as the representative genome 
based on its high completeness (91.64%) and low-contamination (1.36%).
Using Anvio v8, gene annotations were added to based on the Clusters of
Orthologous Groups (COG) database. Functional annotations were used to
extract sequences matching genes of interest, including RecA Recombinase (recA), 
RNA polymerase subunit beta (rpoB), and DNA gyrase subunit B (gyrB), for
downstream comparative and evolutionary analyses. 

Single-nucleotide polymorphism (SNP) and single-codon variant (SCV) 
frequencies were calculated using Anvi’o’s variant detection pipeline. 
Aligned reads were mapped to the reference genome using Bowtie2, and variants 
were identified at each position across all samples. SNP frequencies reflect 
the proportion of reads supporting nucleotide-level variants, while SCV frequencies 
quantify codon-level variation within protein-coding genes. These frequencies were 
extracted per gene and used for statistical analysis of 
site- and time-dependent variation.

## Curation Rational:

This dataset was created to analyze variation in SNPs and SCVs in the
thermophilic bacteria Thermoflexus Hugenholtzii of 3 commonly studied
core genes rpoB- the beta subunit of RNA polymerase in bacteria, gyrB-
the beta subunit of DNA gyrase, and recA- a key protein involved in the
homologous DNA repair process.

    # Datasets
    snp <- read.csv("Thermo_Hugo_SNP.csv")
    scv <- read.csv("Thermo_Hugo_SCV.csv")
    head(scv)

    ##   entry_id unique_pos_identifier                         contig_name     sample_id Temperature    Month Site
    ## 1      398                   108 NODE_120_length_40467_cov_53.543128 June_A_thermo        82.8     June    A
    ## 2      399                   108 NODE_120_length_40467_cov_53.543128 June_C_thermo        73.0     June    C
    ## 3      642                   232 NODE_120_length_40467_cov_53.543128 June_C_thermo        73.0     June    C
    ## 4      106                    23 NODE_120_length_40467_cov_53.543128  Feb_B_thermo        70.0 February    B
    ## 5      107                    23 NODE_120_length_40467_cov_53.543128 June_B_thermo        83.0     June    B
    ## 6      643                   233 NODE_120_length_40467_cov_53.543128 June_C_thermo        73.0     June    C
    ##   ANVIO_Gene_Call Gene_ID   scv_freq
    ## 1             541    rpoB 0.10866373
    ## 2             541    rpoB 0.14428571
    ## 3             541    rpoB 0.08706625
    ## 4             541    rpoB 0.24576271
    ## 5             541    rpoB 0.24277457
    ## 6             541    rpoB 0.08848614

    head(snp)

    ##   entry_id unique_pos_identifier                         contig_name     sample_id Temperature    Month Site
    ## 1      587                   190 NODE_120_length_40467_cov_53.543128 June_A_thermo        82.8     June    A
    ## 2      588                   190 NODE_120_length_40467_cov_53.543128 June_C_thermo        73.0    June     C
    ## 3      753                   284 NODE_120_length_40467_cov_53.543128 June_C_thermo        73.0     June    C
    ## 4      209                    56 NODE_120_length_40467_cov_53.543128  Feb_B_thermo        70.0 February    B
    ## 5      210                    56 NODE_120_length_40467_cov_53.543128 June_B_thermo        83.0     June    B
    ## 6      752                   283 NODE_120_length_40467_cov_53.543128 June_C_thermo        73.0     June    C
    ##   ANVIO_Gene_Call Gene_ID   snp_freq
    ## 1             541    rpoB 0.09397944
    ## 2             541    rpoB 0.13684961
    ## 3             541    rpoB 0.07507886
    ## 4             541    rpoB 0.23430962
    ## 5             541    rpoB 0.22285714
    ## 6             541    rpoB 0.08090957

## Discussion of Biases

Sampling locations in this study were limited and may not capture the
full microbial diversity present in the spring. Thermoflexus
hugenholtzii was selected for analysis due to its high abundance and
consistent presence across all sample sites. The metagenome-assembled
genome (MAG) for T. hugenholtzii was derived from the June Site B
sample, as it had the highest completeness and lowest contamination
among all MAGs generated by MetaBAT2. This MAG was subsequently used as
the reference genome for comparative analyses across all sample sites.

## Personal and Sensitive Information

There is no personal or sensitive information included in these data.
