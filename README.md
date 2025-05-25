# BRCA_immune_classification
This code generates an immune classification based on the xCell matrix derived from RNA-seq data, using a finite-loop immune classification method.
## Method steps
### Step 0
Prepare the xCell data and define the candidate cell types.
### Step 1
Set the number of cluster number (K=2) and the minimal number of select cell types (min_n_cell=6). 
### Step 2
Use the combn function to generate all possible combinations of min_n_cell to length(candidate_celltypes) from the candidate cell types. 
For each combination, the ExecuteCNMF method from the CancerSubtypes package was applied to cluster the samples into two groups.
### Step 3
Evaluated each cluster result by calculating the average silhouette width and the consistency of treatment response within clusters. 
For each cluster result, we assessed the consistency between clusters and treatment response using the average of the maximum proportion of pCR or non-pCR samples within each cluster.
### Final Step
To explore different immune classification models, all possible combinations of 6 to 13 immune cell types were generated from the 13 candidate types, and steps 2 and 3 were repeated for each combination. Only clustering solutions with an average silhouette score > 0.7 and response consistency > 0.6 were considered acceptable. From these, select the combination of cell types that achieves the highest average silhouette width and response consistency as the optimal immune classification model.
## Reference
[1] Xu, T. et al., CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation and visualization. Bioinformatics, 2017. 33(19): p. 3131–3133.
