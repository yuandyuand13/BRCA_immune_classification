My_CNMF_cluster_num = function(n_cell,candidate_celltypes,clusterK=2,xcell_matrix,clinical_data, output_dir){
  # clusterK (int) : the cluster numbers
  # n_cell (int): select n_cell cell types from the candidate cell types
  # candidate_celltypes (character vector) : character vector of candidate cell types
  # xcell_matrix (matrix) : the matrix of xcell result, rows are the cell types and columns are samples
  # clinical_data (data.frame): treatment response of the samples, col.names are c('SampleName','Reponse'), row.names is SampleName
  # output_dir : the output path
  CNMF_seed = 99
  set.seed(CNMF_seed)
  cell_num = paste0(n_cell,'celltypes')
  combinations <- combn(candidate_celltypes, n_cell)
  print(paste(n_cell,"cell types selected from",length(candidate_celltypes),"candidate cell types generating:", ncol(combinations),"combinations."))
  cluster_path3 = file.path(output_dir,paste0(n_cell,'celltypes_',ncol(combinations)))
  dir.create(cluster_path3, recursive = TRUE)
  
  results_tmp = data.frame(matrix(data=0,ncol=5,nrow=0))
  colnames(results_tmp) = c('nCells','m','celltypes','mean_Response','Average_silhouette')
  tmp_row = 1
  m = 1
  for (m in 1:ncol(combinations)){
    num = paste0(n_cell,'celltypes_',m)
    print(num)
    select_celltypes = c(combinations[,m])
    data3 = xcell_matrix[select_celltypes,]
    result1 = tryCatch({
      ExecuteCNMF(data3, clusterNum=clusterK, nrun=30)
    },error=function(e){
      NULL
    })
    
    if (!is.null(result1)) {
      
      print(paste0("n_cell=",n_cell,"_m=",m,": NMF success"))

      group = result1$group
      distanceMatrix = result1$distanceMatrix
      Cluster=as.data.frame(group)
      rownames(Cluster)=colnames(xcell_matrix)
      colnames(Cluster)="Cluster"
      Cluster$Cluster=paste0("C", Cluster$Cluster)
      
      cluster_response = cli_sub[colnames(xcell_matrix),]
      cluster_response$Cluster = factor(Cluster$Cluster,levels = sort(unique(Cluster$Cluster)))
      
      res=My_calculate_Cluster_AverageSilhouette_ResponseConsistency(cluster_response, distanceMatrix)
      tmp = list(n_cell,m,list(select_celltypes),res$mean_Response,res$Average_silhouette)
      results_tmp[tmp_row,] = tmp
      tmp_row = tmp_row+1
      
      if(res$mean_Response>0.6 & res$Average_silhouette>0.7){
        clusterOut=rbind(ID=colnames(Cluster), Cluster)
        write.table(clusterOut, file=file.path(cluster_path3,paste0(num,"_CNMF_seed=",CNMF_seed,".txt")), sep="\t", quote=F, col.names=F)
        
        #对分型的结果进行可视化
        pdf(file=file.path(cluster_path3,paste0(num,"_CNMF_seed=",CNMF_seed,".pdf")), width=10,height=6,onefile = FALSE)
        res=My_plot_Cluster_Silhouette_Response_heatmap(cluster_response, distanceMatrix, similarity=TRUE)
        dev.off()
        save(distanceMatrix,file=file.path(cluster_path3,paste0(num,"_CNMF_seed=",CNMF_seed,".rdata")))
      }
    } else {
      print(paste0("n_cell=",n_cell,"_m=",m,": NMF ERROR"))
      
    }
  }
  print("--------------")
  save(results_tmp,file=file.path(cluster_path3,'results_tmp.rdata'))
  return(results_tmp)
}