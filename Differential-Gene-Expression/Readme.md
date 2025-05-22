# Entropy-Based Differential Gene Expression (DGE)
This folder contains the complete set of scripts used for performing Differential Gene Expression (DGE) analysis using an information-theoretic approach based on Weighted Relative Entropy (WRE). Unlike traditional methods, this approach does not require predefined control groups, making it especially suitable for complex datasets such as antimicrobial resistance (AMR) profiles.

## Method Overview
The method calculates how much the probability distribution of gene expression in each profile deviates from a uniform baseline. It supports:
- General DGE analysis across all profiles
- Test vs control comparisons 
- Output of log2 fold change and statistical significance using Welch’s t-test

## Prerequisites
Before running the scripts, ensure the following R libraries are installed:
- library(tidyverse)
- library(stringr)

## Contents
- human_sample_diff_genes_with_logfc_pval.r – Script for DGE analysis on Human 54 tissue data file. ("human_54samples_TMM_CPM.csv") This produces two CSV files containing overexpressed and underexpressed genes along with their Log2FC and P-value repectively. This scripts performs DGE for TEST vs CONTROL.

- over_and_under_expressed_genenames_fi_by_bi.r – Script which Identifies and exports list of over/underexpressng genes into two separate txt files. This script performs genralized DGE where Test vs Control is not clear.

- plot_fi_by_bi_with_cutoff.r – Script which exports plots for relative expression probability (fi/bi) for profiles for given genes.

- WRE.r – Script which calculated weighted relative entropy for each Gene and exports the WRE values in a CSV file.

