#' KKL clustering
#'
#' @description
#'
#' This function is a wrapper that executes all steps of
#' KKL clustering analysis in one go.
#'
#' @param dataMatrix row is cell, col is feature
#' @param outlier_q the proportion of outliers that need to be removed
#' @param down_n downsampling number
#' @param knn_range the range of the number of neighbors in the KNN graph structure
#' @param iter the number of iterations
#' @param compute_index the value of clustering quality index that need to be calculated
#' @param assess_index evaluation index used to select the optimal KNN graph structure
#' @param cores  the number of threads
#' @param seed random seed
#'
#' @export
#'
kkl <- function(dataMatrix = NULL,
                outlier_q = 0.2,
                down_n = 300,
                knn_range = 5:70,
                iter = 50,
                compute_index =  c("Davies_Bouldin","Calinski_Harabasz"),
                assess_index = "Davies_Bouldin",
                cores = 1,
                seed = 723){

  outlier_kmeans = OutlierKmeans(dataMatrix = dataMatrix,
                                 outlier_q = outlier_q,
                                 down_n = down_n,
                                 cores = cores,
                                 seed = seed)

  # random sampling and choose optmial KKN graph structure
  sampling_result = SamplingLouvain(dataMatrix = dataMatrix,
                                    outlier_kmeans = outlier_kmeans,
                                    knn_range = knn_range,
                                    iter = iter,
                                    compute_index = compute_index,
                                    cores = cores,
                                    seed = seed)

  # get final clustering result
  new_louvain = NewLouvain(sampling_result = sampling_result,
                           outlier_kmeans = outlier_kmeans,
                           assess_index = assess_index,
                           cores = cores)

  cluster_df = data.frame(row.names = rownames(new_louvain$cluster_meta),
                          name = rownames(new_louvain$cluster_meta),
                          cluster = as.factor(new_louvain$cluster_meta$cluster),
                          stringsAsFactors = FALSE)

  result = list(outlier_kmeans = outlier_kmeans,
                sampling_result = sampling_result,
                new_louvain = new_louvain,
                cluster_df = cluster_df)

  return(result)

}
