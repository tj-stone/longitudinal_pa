# Longitudinal pilocytic astrocytoma R scripts

R scripts for processing of Illumina methylation array data and RNA-seq data plus generation of associated plots. Used in "The tumour microenvironment of pilocytic astrocytoma evolves over time via enrichment for microglia." by Stone et al.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Hardware requirements](#hardware-requirements)
- [Software Requirements](#software-requirements)

# Overview

`longitudinal_PA_RNA_Seq.R` contains the step by step commands used to read in and process RNA-Seq data, plus subsequent quality control filtering, differential expression analysis, and gene set variation analysis for enrichment against cell type-specific gene sets.

`longitudinal_PA_methylation.R` contains the step by step commands used to read in Illumina EPIC methylation array data, perform quality control and pre-processing, plus subsequent differential methylation analysis and over-representation analysis against KEGG, MSigDB Hallmark, and cell type-specific gene sets.

`longitudinal_PA_oncoprint.R` contains the code used to generate a summary oncoprint plot of mutation data for variants identified within the cohort.

# System Requirements
## Hardware Requirements
These scripts should be compatible with any standard system that supports *R 4.1.1* with enough RAM to support the constituent operations. Scripts were implemented originally on a system with 8Gb RAM.

## Software Requirements
The scripts provided depend on a number of R packages (listed within each script). Functionality for all scripts has been tested on *R 4.1.1* operating under *macOS Catalina 10.15.7*
