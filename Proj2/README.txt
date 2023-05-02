This is the project doing the second data analysis on MS data.

We have the data about MS of in vitro T cells for SNS101.

The data are in ../Data

For this first try, we do heatmap and pca.

#one issue related to do normalize by quantile method,"return code from pthread_create() is 22"

To fix: BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)
and then need to restart!!!!

#Note:
Curve fitting is largely dependent on the number of consecutive measurements. In proteomics, missing values
are not uncommon. If too many points are missing, a curve cannot be fit accurately. This is especially
common for small batches. In this case, we suggest to not Åt the curve to the speciÅc peptide within the
speciÅc batch, and proceed directly to discrete correction methods. To identify such peptides, absolute and
relative thresholds (abs_threshold and pct_threshold) on the number of missing values for each peptide
can be passed as parameters to adjust_batch_trend_df().
