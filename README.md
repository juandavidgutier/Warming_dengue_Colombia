# Project Title

Warming-sensitive dengue risk emerges at higher elevations in Colombia: evidence from Stochastic Treatment Regimes

## Description

Code and dataset shared to reproduce the results of the paper Warming-sensitive dengue risk emerges at higher elevations in Colombia: evidence from Stochastic Treatment Regimes.
The file data\_final.csv is the dataset used for the results presented in the manuscript.
To obtain the adjustment set required to control for confounding use the file DAG.py. To reproduce the results of the STR-TMLE framework use the file Str\_tmle.r. To reproduce the results of the g-computation use the file logistic\_gcomputation.R

## Data Privacy and Anonymization

This dataset has been processed to ensure complete anonymization and contains no personally identifiable information (PII). All data has been:

* Aggregated at appropriate spatial/temporal scales
* Stripped of any individual identifiers
* Processed to remove direct or indirect identifying elements

The dataset is suitable for public sharing and complies with data privacy standards.

## Privacy Statement

This repository contains datasets that have been carefully processed to protect individual privacy:

### What is NOT included:

* Names, addresses, or contact information
* Individual-level identifiers
* Location data below 25 km resolution
* Timestamps more precise than monthly
* Any data that could be used to re-identify individuals

### Data Processing:

* Spatial aggregation to municipality
* Temporal aggregation to monthly averages

### Compliance:

This dataset meets requirements for public data sharing under applicable privacy regulations.

## Author

Juan David Gutiérrez

## libraries

haldensify; sl3; tmle3; tmle3shift; dplyr; ggplot2; caret
pandas; numpy; dowhy; statsmodels; matplolib; scipy;

