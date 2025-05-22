# Jensen-Shannon Distance & Mutual Information Analysis
This folder contains scripts for computing Jensen-Shannon (JS) Distance, two-variable mutual information (MI), and three-variable multivariate mutual information (MMI) across gene expression profiles. These information-theoretic measures are used to uncover both pairwise and higher-order relationships between antibiotic resistance phenotypes.

## Method Overview
### Jensen-Shannon (JS) Distance:
Measures the divergence between probability distributions of gene expression across profiles.

### Two-variable Mutual Information (MI):
Captures non-linear dependencies between two profiles, helping reveal co-expression or shared resistance mechanisms.

### Three-variable Mutual Information (MMI):
Identifies profiles with redundant or synergistic information, useful for feature selection and interdependency analysis.

## Prerequisites
Before running the scripts, ensure the following R libraries are installed:
- library(tidyverse)
- library(stringr)
- library(ggplot2)
- library(reshape2)
- library(scales)


## Contents
- js_distances_between_samples.r – Computes pairwise JS distance across all resistance profiles for "Japan_Lab_Transcriptomic_Data1.csv" and export a csv file along with a heatmap.

- mutual_info_data1.r – Calculates pairwise MI across all resistance profiles for "Japan_Lab_Transcriptomic_Data1.csv" and export a csv file along with a heatmap.

- mutual_info_data2.r – Calculates pairwise MI across all resistance profiles for "Japan_Lab_Transcriptomic_Data2.csv" and export a csv file along with a heatmap.

- mmi_data1.r – Performs three-variable mutual information analysis for "Japan_Lab_Transcriptomic_Data1.csv" and export a csv file along with  optional barplot.

- mmi_data2.r – Performs three-variable mutual information analysis for "Japan_Lab_Transcriptomic_Data2.csv" and export a csv file along with  optional barplot.
