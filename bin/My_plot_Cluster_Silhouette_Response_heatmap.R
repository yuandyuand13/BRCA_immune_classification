My_plot_Cluster_Silhouette_Response_heatmap = function (cluster_response, distanceMatrix = NULL, similarity = TRUE) 
{
  if(!all(colnames(distanceMatrix) == rownames(cluster_response))){
    stop(paste0("distanceMatrix samples and cluster samples do not match"))
  }
  mean_Response = mean(apply(prop.table(table(cluster_response$Cluster,cluster_response$Response),margin = 1), 1, max))
  
  group = as.numeric(gsub("[a-zA-Z]",'',cluster_response$Cluster))
  clusterNum = length(unique(group))
  if (!is.null(distanceMatrix[1, 1])) {
    layout(matrix(c(1, 2), 1, 2, byrow = FALSE), widths = c(6, 4))
  }
  if (class(distanceMatrix) == "Similarity") {
    si = silhouette_SimilarityMatrix(group, distanceMatrix)
  } else {
    si = silhouette(group, distanceMatrix)
  }
  Silhouette = as.data.frame.array(si)
  Average_silhouette = mean(Silhouette$sil_width)  
  
  
  ind = order(group, -si[, "sil_width"])
  
  cluster = cluster_response[ind,c("Cluster","Response")]

  library(stringr)
  library(circlize)
  library(RColorBrewer)
  hhh = cluster
  col_map = list(clusterCol = c("#D53A58","#54CB3E","#1F83E0","#7700FF","#FF9900"),
                 bioCol2 = brewer.pal(5, 'Set2'),
                 bioCol3 = brewer.pal(5, 'Set3'),
                 bioCol4 = brewer.pal(5, 'Set1'),
                 bioCol5 = brewer.pal(5, 'Paired'),
                 bioCol5 = brewer.pal(5, 'Accent'))
  anno_color = list()
  for(i in 1:ncol(hhh)){
    col_name = colnames(hhh)[i]
    Col = col_map[[i]]
    col_color = Col[1:length(unique(hhh[,i]))]
    names(col_color) = sort(unique(hhh[,i]))
    anno_color[[i]] = col_color
  }
  names(anno_color) = colnames(hhh)
  
  attr(distanceMatrix, "class") = NULL
  consensusmap(distanceMatrix, 
               Rowv = ind, Colv = ind,
               main = "Clustering display",
               annCol = cluster_response[,c("Cluster","Response")], 
               annColors = anno_color, 
               labRow = NA, 
               labCol = NA, 
               cexRow = 2,
               scale = "none")
  plot(si, col = anno_color$Cluster,oma = c(0,0, 0, 0))
  
  
  par(mfrow = c(1, 1))
  
  out <- list(mean_Response=mean_Response, Average_silhouette=Average_silhouette)
  return(out)
}
