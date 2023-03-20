# Microbiome analysis

This is the home for the bioinformatic and statistical analysis of the microbiome analysis.

## 1_start

This is where meta data are saved.

## 2_data

This is where sequencing raw data were stored. Before running the bioinformatic pipeline, place your data here. The runs must be named *runXXXX_r1.fastq.gz* and *runXXXX_r2.fastq.gz* where *XXXX* correspond to the run-number from the design-file. The data can be found on european nucleotide archive under the accession **PRJEB59165** .

## 3_pipeline

This is the bioinformatic pipeline. Calculations were performed at sciCORE (<http://scicore.unibas.ch/>) scientific computing center at University of Basel.

**01: Demultiplexing:** This folder contains the code for demultiplexing raw data.\

**02: ASV dadapipe:** This folder contains all code to run the dadapipe. With dada2 sequences were quality filtered, truncated, dereplicated and denoised. A count table was created and taxonomy was assigned. DECIPHER was used to cluster sequences by similarity.

## 4_output

This is where the pipeline output is stored.

## 5_stats

This is the statistical analysis of the microbiome. The folder \*\*2_scripts contains all used scripts:

**Analysis scripts:** The script Microbiomes1_data_preparation.Rmd imports microbiome data and normalize them. Microbiomes2_maize_harvest.Rmd analyses data after the conditioning phase, Microbiomes3_wheat_sowing.Rmd during sowing phase and Microbiomes4_wheat_growth.Rmd at the end of the wheat growth phase.\

**Microbiomes0_wrapper.Rmd:** By knitting this script you receive one pdf out of all child Rmd-documents.
