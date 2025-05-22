# Analysis of Antimicrobial Resistance Data Using Information Theory

This repository contains the complete workflow, analysis scripts, and models for studying Antimicrobial Resistance (AMR) using information-theoretic methods. The project integrates entropy-based gene expression analysis, feature selection through mutual information, and machine learning to enable robust prediction of resistance profiles.
Developed as part of an academic research project, this framework is designed to be both scalable and interpretable, especially in settings where classical statistical assumptions—like defined control groups—do not hold.


## 📁 Repository Structure
- .
- ├── Differential-Gene-Expression/                # Entropy-based DGE using Weighted Relative Entropy
- ├── Js-Distance-&-Mutual-Information/            # JS distance, pairwise and multivariate mutual information
- ├── AMR-Prediction-Using-XGBoost/                # AMR prediction using XGBoost with different feature sets
- ├── SuperInformation/                            # Higher-order entropy analysis of genomic sequences


## Key Features
### Entropy-Based Differential Gene Expression (DGE):
Identifies differentially expressed genes using Weighted Relative Entropy (WRE)—works with or without defined test/control groups.

### Mutual Information & Jensen-Shannon Distance:
Quantifies non-linear dependencies and distributional differences between resistance profiles.

### Multivariate Mutual Information (MMI):
Selects non-redundant, high-information features for efficient AMR prediction.

### Machine Learning with XGBoost:
Predicts AMR expression profiles using minimal feature sets with performance close to full models.

### Superinformation (Experimental):
Applies higher-order entropy to genomic sequences to identify complex resistance patterns beyond pairwise relationships.


## 🧰 Requirements
### R Packages
Used in entropy-based DGE and mutual information scripts:
- library(tidyverse)
- library(stringr)
- library(ggplot2)
- library(scaler)
- library(reshape2)

### Python Packages
Used for AMR prediction via machine learning:
- import pandas as pd
- import numpy as np
- import matplotlib.pyplot as plt
- import seaborn as sns
- from xgboost import DMatrix, train
- from sklearn.model_selection import train_test_split
- from sklearn.metrics import mean_squared_error


## 📂 Folder Highlights

### Differential-Gene-Expression/
Implements Weighted Relative Entropy for differential gene expression. Works for both general and test/control setups.

### Js-Distance-&-Mutual-Information/
Scripts to calculate:
- JS Distance for profile divergence
- Mutual Information (MI) for pairwise gene correlation
- Multivariate MI (MMI) for feature selection and profile uniqueness

### AMR-Prediction-Using-XGBoost/
Collab notebooks comparing:
- All features
- MMI-selected principal features
- All features excluding principal profiles
- Includes comparision_graph.ipynb to visualize R² scores.

### SuperInformation/
Performs exploratory analysis on DNA sequences using superinformation—a entropy measure that detects higher-order dependencies.


## Acknowledgements
**This research was conducted under the guidance of Dr. Nithya Ramakrishnan and supported by the Department of IT, BT and S&T, Govt. of Karnataka. Thanks to the lab team for their input and support.**
