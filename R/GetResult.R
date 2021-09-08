#' Louvain clustering with optimal KNN structure
#'
#' @description
#'
#' This function is used to run Louvain clustering with optimal KNN graph structure
#'
#' @param sampling_result sampling result from SamplingLouvain function
#' @param outlier_kmeans outlier detection result from OutlierKmeans function
#' @param assess_index evaluation index used to select the optimal KNN graph structure
#' @param cores  the number of threads
#'
#' @export
#'
NewLouvain <- function(sampling_result = NULL,
                       outlier_kmeans = NULL,
                       assess_index = "Davies_Bouldin",
                       cores = 1){

  index = assess_index

  int_best = sampling_result$int_best
  ds_kmeans = outlier_kmeans$ds_kmeans
  ds_kmeans_centers = outlier_kmeans$ds_kmeans_new_centers
  knn_best = as.numeric(str_split_fixed(int_best[index,"best"],"k",2)[,2])

  cat("\n/// louvain with", knn_best, "\n")
  d_mick = parallelDist::parDist(x = ds_kmeans_centers, method = "euclidean", diag = TRUE, upper = TRUE, threads = 10)
  knn = dbscan::kNN(d_mick, k = knn_best)
  knn_edge = reshape2::melt(knn$id)
  knn_graph = graph_from_edgelist(as.matrix(knn_edge[,c(1,3)]), directed = FALSE)
  c_louvain = cluster_louvain(knn_graph)

  cat("\n/// mapping \n")
  c_louvain_cluster = c_louvain$membership
  names(c_louvain_cluster) = 1:length(c_louvain_cluster)
  gc()
  map_cluster = GetMapClusterLabel(orig_cluster = ds_kmeans$cluster,
                                   second_cluster = c_louvain_cluster,
                                   cores = cores)
  cluster_meta = data.frame(row.names = names(map_cluster),
                            cluster = map_cluster,
                            stringsAsFactors = FALSE)

  result = list(knn_best = knn_best,
                louvain_result = c_louvain,
                cluster_meta = cluster_meta)

  return(result)
}

#' Map the K-means results to all cells
#'
#' @Description:
#' This function is a built-in function of NewLouvain,
#' which used to map clustering results to all cells
#'
#' @param orig_cluster clustering label of first time
#' @param second_cluster clustering label of second time
#' @param cores  the number of threads
#'
#' @export
#'
GetMapClusterLabel <- function(orig_cluster = NULL,
                               second_cluster = NULL,
                               cores = 1){
  orig = orig_cluster
  new = second_cluster
  max_cores = detectCores()-1

  if(cores==1 | cores<max_cores){
    tmp = lapply(1:length(new), function(i){
      n = names(new)[i]
      ind = which(as.character(orig)==n)
      tmp = rep(new[n],length(ind))
      names(tmp) = names(orig[ind])
      tmp })
  }else{
    cl = makeCluster(cores)
    clusterExport(cl = cl, varlist = c("orig","new"), envir=environment())
    tmp = parLapply(cl, 1:length(new),
                    function(i){
                      n = names(new)[i]
                      ind = which(as.character(orig)==n)
                      tmp = rep(new[n],length(ind))
                      names(tmp) = names(orig[ind])
                      tmp
                    })
    stopCluster(cl)
  }

  final = unlist(tmp)[names(orig)]
  return(final)
}
