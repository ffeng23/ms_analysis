#R code for testing proBatch code

require(dplyr)
require(tibble)
require(ggplot2)


feature_id_col = 'peptide_group_label'
measure_col = 'Intensity'
sample_id_col = 'FullRunName'
essential_columns = c(feature_id_col, measure_col, sample_id_col)

technical_factors = c('MS_batch', 'digestion_batch')
biological_factors = c('Strain', 'Diet', 'Sex')
biospecimen_id_col = 'EarTag'

batch_col = 'MS_batch'

library(proBatch)
data('example_proteome', 'example_sample_annotation', 'example_peptide_annotation',
package = 'proBatch')

generated_sample_annotation <- date_to_sample_order(example_sample_annotation,
time_column = c('RunDate','RunTime'),
new_time_column = 'generated_DateTime',
dateTimeFormat = c('%b_%d', '%H:%M:%S'),
new_order_col = 'generated_order',
instrument_col = NULL)

library(knitr)
kable(generated_sample_annotation[1:5,] %>%
select(c('RunDate', 'RunTime', 'order', 'generated_DateTime', 'generated_order')))

generated_peptide_annotation <- create_peptide_annotation(example_proteome,
feature_id_col = 'peptide_group_label',
protein_col = 'Protein')
example_proteome = example_proteome %>% select(one_of(essential_columns))
gc()


example_matrix <- long_to_matrix(example_proteome,
feature_id_col = 'peptide_group_label',
measure_col = 'Intensity',
sample_id_col = 'FullRunName')

log_transformed_matrix <- log_transform_dm(example_matrix,
log_base = 2, offset = 1)

color_list <- sample_annotation_to_colors(example_sample_annotation,
factor_columns = c('MS_batch', 'digestion_batch',
'EarTag', 'Strain',
'Diet', 'Sex'),
numeric_columns = c('DateTime','order'))


plot_sample_mean(log_transformed_matrix, example_sample_annotation, order_col = 'order',
batch_col = batch_col, color_by_batch = TRUE, ylimits = c(12, 16.5),
color_scheme = color_list[[batch_col]])

log_transformed_long <- matrix_to_long(log_transformed_matrix)
batch_col = 'MS_batch'
plot_boxplot(log_transformed_long, example_sample_annotation,
batch_col = batch_col, color_scheme = color_list[[batch_col]])

median_normalized_matrix = normalize_data_dm(log_transformed_matrix,
normalize_func = 'medianCentering')

quantile_normalized_matrix = normalize_data_dm(log_transformed_matrix,
normalize_func = 'quantile')

plot_sample_mean(quantile_normalized_matrix, example_sample_annotation,
color_by_batch = TRUE, ylimits = c(12, 16),
color_scheme = color_list[[batch_col]])

selected_annotations <- c('MS_batch',
'digestion_batch', 'Diet')

plot_hierarchical_clustering(quantile_normalized_matrix,
sample_annotation = example_sample_annotation,
color_list = color_list,
factors_to_plot = selected_annotations,
distance = 'euclidean', agglomeration = 'complete',
label_samples = FALSE)

plot_heatmap_diagnostic(quantile_normalized_matrix, example_sample_annotation,
factors_to_plot = selected_annotations,
cluster_cols = TRUE,
color_list = color_list,
show_rownames = FALSE, show_colnames = FALSE)

pca1 = plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = 'MS_batch',
plot_title = 'MS batch', color_scheme = color_list[['MS_batch']])
pca2 = plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = 'digestion_batch',
plot_title = 'Digestion batch', color_scheme = color_list[['digestion_batch']])
pca3 = plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = 'Diet',
plot_title = 'Diet', color_scheme = color_list[['Diet']])
pca4 = plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = 'DateTime',
plot_title = 'DateTime', color_scheme = color_list[['DateTime']])

library(ggpubr)
ggarrange(pca1, pca2, pca3, pca4, ncol = 2, nrow = 2)

pca_spec = plot_PCA(quantile_normalized_matrix, example_sample_annotation,
color_by = 'digestion_batch',
plot_title = 'Digestion batch')

pca_spec

#this might take very long time
plot_PVCA(quantile_normalized_matrix, example_sample_annotation,
technical_factors = technical_factors,
biological_factors = biological_factors)

quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix)

loess_fit_df <- adjust_batch_trend_df(quantile_normalized_long, example_sample_annotation)

loess_fit_30 <- adjust_batch_trend_df(quantile_normalized_long, example_sample_annotation,
span = 0.3)
plot_with_fitting_curve(feature_name = '10231_QDVDVWLWQQEGSSK_2',
fit_df = loess_fit_30, fit_value_col = 'fit',
df_long = quantile_normalized_long,
sample_annotation = example_sample_annotation,
color_by_batch = TRUE, color_scheme = color_list[[batch_col]],
plot_title = 'Span = 30%')

loess_fit_70 <- adjust_batch_trend_df(quantile_normalized_long, example_sample_annotation,
span = 0.7)
plot_with_fitting_curve(feature_name = '10231_QDVDVWLWQQEGSSK_2',
fit_df = loess_fit_70, fit_value_col = 'fit',
df_long = quantile_normalized_long,
sample_annotation = example_sample_annotation,color_by_batch = TRUE, color_scheme = color_list[[batch_col]],
plot_title = 'Span = 70%')


peptide_median_df <- center_feature_batch_medians_df(loess_fit_df, example_sample_annotation)
plot_single_feature(feature_name = '10231_QDVDVWLWQQEGSSK_2', df_long = peptide_median_df,
sample_annotation = example_sample_annotation, measure_col = 'Intensity',
plot_title = 'Feature-level Median Centered')

comBat_df <- correct_with_ComBat_df(loess_fit_df, example_sample_annotation)

plot_single_feature(feature_name = '10231_QDVDVWLWQQEGSSK_2',
df_long = loess_fit_df,
sample_annotation = example_sample_annotation,
plot_title = 'Loess Fitted')
plot_single_feature(feature_name = '10231_QDVDVWLWQQEGSSK_2',
df_long = comBat_df,
sample_annotation = example_sample_annotation,
plot_title = 'ComBat corrected')