#### Description

- boruta: Data which dimensions were reduced by selecting most relevant descriptors/attributes  based on the feature selection algorithm Boruta (R).
- cfs: Data which dimensions were reduced by selecting most relevant descriptors/attributes  based on the feature selection algorithm CFS (Correlation-based Feature Selection) algorithms (Weka).
- mean: Data which dimensions were reduced by aggregating each property's 10 associated descriptors using the mean. Thus, there are only 48 attributes, 1 per property.
- pca: Data which dimensions were reduced by using only the 30 first principal components (>95% explained variance) of PCA applied to original data.
- pca_per_property: Data which dimensions were reduced by aggregating each property's 10 associated descriptors using PCA. Thus, data has 48 attributes, 1 per property. The attribute for the property k0 is the first principal component of PCA applied to AASA1K0, ..., AASA10K0.
