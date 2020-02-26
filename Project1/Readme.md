## Comparative study on PD-L1 and PD-1 in cancers

This is a #dailycoding project. I utilized the "FireBrowse" resource (via firebrowseR package in R) to access the
expression data across all TCGA cohorts. I obtained the median log2 expression level of PD-L1 and PD-1 and
compared using a scatter plot (using ggplot2 package) across TCGA cohorts. I observed that these two genes
show a significant linear relationship across all cohorts (p=0.002, r2 = 0.24). It also suggests their consequential
role in promoting cancer and based on this analysis, DLBC, PAAD, and BRCA could be considered as the most
relevant oncological disease cohorts for their solid tumor immunotherapy.
