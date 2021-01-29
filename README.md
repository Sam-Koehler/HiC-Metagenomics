# Metagenomic pipeline using Shotgun and Hi-C reads

This pipeline is based in Snakemake and analyzes a multi-sample metagenomic dataset of FASTQ reads, ultimately outputting graphs of both the taxonomic community profile and percentage abundances of species present.


# Motivation

Metagenomics poses unique and challenging problems in comparison to the standard study of single organism datasets, for which most sequencing platforms are designed. Hi-C proximity ligation allows for improved deconvolution of metagenomic sequencing data by providing information about the co-localization of sequences within a single cell. 


# Data Overview

The data used in this pipeline consists of four samples of FASTQ read data of the algal species *Nannochloropsis oceanica*. *N. oceanica* is of significant interest as an alternative energy source due to its high lipid content and minimal growth requirements. However, *N. oceanica* is easily susceptible to infection by microscopic pathogen(s) that are not yet identified. By using metagenomic techniques and leveraging the increased deconvolutional power of Hi-C, we can identify the organisms present in algal ponds. We are studying four samples, two of which represent healthy ponds, and two that represent sick ponds that are no longer viable for energy production. In this project we aim to identify the species present in each of the ponds, and upon comparison of the healthy and sick ponds, identify the possible pathogen(s).


# Future Directions

In order to make this pipleine usable for other metagenomic datasets, we aim to create an open-source UI that allows for input of the number of samples, the conditions of each sample, and the types of read data. The user will be able to recreate a similar metagenomic analysis by downloading the Snakemake pipeline and the required environments.


# Acknowledgements

All FASTQ data is provided by Los Alamos National Lab through the University of Arizona. This project was created at the University of Oregon. 
