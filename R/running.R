#' CDSKNN clustering
#'
#' @description
#'
#' This function is a wrapper that executes all steps of
#' CDSKNN clustering analysis in one go.
#'
#' @param dat expression matrix, row is cell, col is feature
#' @param partition_count the number of partitions
#' @param batch_size the number of threads
#' @param outlier_det whether running outlier detection
#' @param outlier_methods outlier detection methods, md or ed.
#' "md" utilizes the Mahalanobis distance for detection,
#' while "ed" treats points that are far from the center proportionally as outliers.
#' @param outlier_q the proportion of outliers you want to filter.
#' When outlier_methods is ed, this parameter takes effect.
#' @param min_n the minimum number of data points for outlier detection.
#' @param num_init the proportion of outliers that need to be removed
#' @param max_iters the maximum iterations for kmeans
#'
#' @param knn_range the range of the number of neighbors in the KNN graph structure
#' @param cluster_method louvain or leiden in resampling
#' @param resolution clustering parameters settings in resampling
#' @param iter the number of iterations in resampling
#' @param is_weight Whether to use distance weights when constructing the KNN graph
#'
#' @param assess_index evaluation index used to select the optimal KNN graph structure
#' @param new_cluster_method louvain or leiden
#' @param res_range resolution range in clustering

#' @param seed random seed
#'
#' @export
#'
CDSKNN <- function(dat = NULL,

                   partition_count = 300,
                   batch_size=500,

                   outlier_det = TRUE,
                   outlier_methods = "md",
                   outlier_q = 0.2,
                   min_n = 200,
                   num_init = 10,
                   max_iters = 1000,

                   cluster_method = "louvain",
                   resolution = 1,
                   knn_range = c(3:10),
                   iter = 20,
                   is_weight = TRUE,

                   assess_index = "Calinski_Harabasz",
                   res_range = seq(0.2,3,0.2),
                   new_cluster_method = "louvain",

                   seed = 723){

  # data("pca_result")
  # dat = pca_result
  # partition_count = 300
  # batch_size=500
  #
  # outlier_det = TRUE
  # outlier_methods = "md"
  # outlier_q = 0.2
  # min_n = 200
  #
  # num_init = 10
  # max_iters = 1000
  #
  # cluster_method = "louvain"
  # resolution = 1
  # knn_range = c(3:50)
  # iter = 50
  # is_weight = TRUE
  #
  # assess_index = "Calinski_Harabasz"
  # res_range = seq(0.2,3,0.2)
  # new_cluster_method = "louvain"
  #
  # seed = 723
  # cores = 1

  cat("\n\n ///// 1) partitioning with outlier detection \n\n")

  outlier_kmeans = pre_partitioning(dat = dat,

                                    partition_count = partition_count,
                                    batch_size = batch_size,

                                    outlier_det = outlier_det,
                                    outlier_methods = outlier_methods,
                                    outlier_q = outlier_q,
                                    min_n = min_n,

                                    num_init = num_init,
                                    max_iters = max_iters,

                                    seed = seed)


  cat("\n\n ///// 2) Finding the optimal graph structure \n\n")

  sampling_result = resampling(dat = dat,
                               outlier_kmeans = outlier_kmeans,
                               resolution = resolution,
                               cluster_method = cluster_method,
                               knn_range = knn_range,
                               iter = iter,
                               is_weight = is_weight,
                               seed = seed)

  cat("\n\n ///// 3) Finding the optimal resolution \n\n")

  new_r = new_clustering(sampling_result = sampling_result,
                         outlier_kmeans = outlier_kmeans,
                         assess_index = assess_index,
                         cluster_method = new_cluster_method,
                         res_range = res_range,
                         is_weight = is_weight,
                         seed = seed)

  result <- list(outlier_kmeans = outlier_kmeans,
                 sampling_result = sampling_result[-1],
                 cluster_result = new_r)

  command = GetCommand()
  command$params = command$params[-1]
  result$command = command
  return(result)

}


#' GetCommand
#'
#' @noRd
#'
GetCommand <- function(){
  time.stamp <- Sys.time()
  command.name = sys.calls()[[1]]
  command.name = strsplit(as.character(command.name[[1]]),"\\(")[[1]]

  argnames <- argnames <- names(x = formals(fun = sys.function(which = sys.parent(n = 1))))
  params <- list()
  p.env <- parent.frame(n = 1)
  argnames <- intersect(x = argnames, y = ls(name = p.env))
  # argnames <- setdiff(argnames, c("mamba","murp","metadata","object"))
  for (arg in argnames) {
    param_value <- get(x = arg, envir = p.env)
    # if (inherits(x = param_value)) {
    #   next
    # }
    params[[arg]] <- param_value
  }

  # p = list(name = command.name, time.stamp = time.stamp, argnames = argnames, params = params)
  p = list(name = command.name, time.stamp = time.stamp, params = params)
  return(p)
}
