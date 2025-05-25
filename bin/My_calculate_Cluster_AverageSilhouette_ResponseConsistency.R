My_calculate_Cluster_AverageSilhouette_ResponseConsistency =function (cluster_response, distanceMatrix = NULL, similarity = TRUE) 
{
  if(!all(colnames(distanceMatrix) == rownames(cluster_response))){
    stop(paste0("distanceMatrix samples and cluster samples do not match"))
  }
  mean_Response = mean(apply(prop.table(table(cluster_response$Cluster,cluster_response$Response),margin = 1), 1, max))
  group = as.numeric(gsub("[a-zA-Z]",'',cluster_response$Cluster))
  if (class(distanceMatrix) == "Similarity") {
    si = silhouette_SimilarityMatrix(group, distanceMatrix)
  } else {
    si = silhouette(group, distanceMatrix)
  }
  Silhouette = as.data.frame.array(si)
  Average_silhouette = mean(Silhouette$sil_width)
  out <- list(mean_Response=mean_Response, Average_silhouette=Average_silhouette)
  return(out)
}
