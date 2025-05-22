# AMR Prediction Using XGBoost
This folder contains Collab notebook codes used to train and evaluate machine learning models for predicting antimicrobial resistance (AMR) profiles using gene expression data. Models are built using the XGBoost regression algorithm, with a focus on comparing different feature selection strategies derived from Multivariate Mutual Information (MMI).

## Method Overview
### Feature Sets Compared:
- Principal Features: Top profiles selected using MMI (informative and non-redundant).
- All Features: Full set of resistance profiles.
- All Features Excluding Principal: All profiles except the selected principal ones.

## Modeling Approach:
- XGBoost regression models are trained using each feature set.
- Performance is evaluated using the R² score and Mean Squared Error (MSE).
- Results are visualized to compare model accuracy across feature sets.

## Prerequisites
Ensure the following Python libraries are installed before running the notebooks:
- import pandas as pd
- import numpy as np
- import matplotlib.pyplot as plt
- import seaborn as sns
- from xgboost import DMatrix, train
- from sklearn.model_selection import train_test_split
- from sklearn.metrics import mean_squared_error

## Contents
- xgboost_principal_features.ipynb – Trains the model using only principal (MMI-selected) features

- xgboost_all_features.ipynb – Trains the model using all resistance profiles

- xgboost_all_features_without_principal_features.ipynb – Trains the model using all resistance profiles but principal (MMI-selected) features are removd from the feature set.

- comparision_graph.ipnyb - Plots the R² values for all the features compbinations in a single graph for better visualization.

