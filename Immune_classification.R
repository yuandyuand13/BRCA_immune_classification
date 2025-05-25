rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd = getwd()
# packcage
library(limma)
library(ggplot2)
library(factoextra)
library(NbClust)
library(CancerSubtypes)
library(NMF)
library(dplyr)
library(xCell)
library(survival)
library(survminer)

source_path = file.path(wd,'bin')
source(file.path(source_path,"My_plot_Cluster_Silhouette_Response_heatmap.R"))
source(file.path(source_path,"My_calculate_Cluster_AverageSilhouette_ResponseConsistency.R"))
source(file.path(source_path,"My_CNMF_cluster_num.R"))


data_path = '../01_data_xcell_process'
load(file=file.path(data_path,'BC.rdata'))
cli = BC[["all"]][["meta_all"]]
xcell = BC[["all"]][["xcell_all"]][1:64,]


candidate_celltypes = c("B-cells", "CD4+ naive T-cells", "CD8+ T-cells", "CD8+ Tcm", 
                        "Class-switched memory B-cells", "ly Endothelial cells", "Mast cells", 
                        "Monocytes", "naive B-cells", "NKT", "Plasma cells", "pro B-cells", 
                        "Tgd cells")

# clinical data 
timeline = 'baseline'
cli_sub = cli[cli$Timeline==timeline,c("SampleName","Response")]
rownames(cli_sub) = cli_sub$SampleName

# xcell matrix
xcell_sub = xcell[candidate_celltypes,rownames(cli_sub)] 

#### start immune classification by finite-loop immune classification method ####
# set the output path
cluster_path = 'Immune_classification_result'
dir.create(cluster_path)
# 01. set the number of clusters (K=2) and the initial minimal number of cells (n_cell=6) ####
K = 2
num_cell_list = c(6:length(candidate_celltypes)) 
num_cell_list = 13
result_all = data.frame()
for (n_cell in num_cell_list) {
  # 02. generate all possible combinations of n_cell cell types from the 13 candidate cell types ####
  # 03. calculating the average silhouette width and the consistency of treatment response within clusters ####
  results_tmp = My_CNMF_cluster_num(n_cell,candidate_celltypes,clusterK=K,
                                    xcell_matrix=xcell_sub,clinical_data=cli_sub,
                                    output_dir=cluster_path)
  result_all = dplyr::bind_rows(result_all, results_tmp)
}

# find the best immune classification result and cell types
results_all_sub = subset(results_all,mean_Response>0.65 & Average_silhouette>0.70)
best_cluster = results_all_sub[results_all_sub$mean_Response==max(results_all_sub$mean_Response),]
select_celltypes = best_cluster$celltypes[[1]]
