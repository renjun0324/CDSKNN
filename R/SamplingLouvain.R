#' Random sampling to do KNN-Louvain clustering
#'
#' @description
#' This function is used to randomly sample each cluster multiple times
#' to select the most powerful KNN graph structure.
#'
#' @param dataMatrix row is cell, col is feature
#' @param outlier_kmeans the result from OutlierKmeans function
#' @param knn_range the range of the number of neighbors in the KNN graph structure
#' @param iter the number of iterations
#' @param compute_index the value of clustering quality index that need to be calculated
#' @param cores the number of threads
#' @param seed random seed
#'
#' @details The currently available indices are:
#' \itemize{
#'    \item{Ball_Hall}\item{Banfeld_Raftery}\item{C_index}\item{Calinski_Harabasz}
#'    \item{Davies_Bouldin}\item{Det_Ratio}\item{Dunn}\item{Gamma}\item{G_plus}
#'    \item{Ksq_DetW}\item{Log_Det_Ratio}\item{Log_SS_Ratio}\item{McClain_Rao}
#'    \item{PBM}\item{Point_Biserial}\item{Ray_Turi}\item{Ratkowsky_Lance}
#'    \item{Scott_Symons}\item{SD_Scat}\item{SD_Dis}\item{S_Dbw}\item{Silhouette}
#'    \item{Tau}\item{Trace_W}\item{Trace_WiB}\item{Wemmert_Gancarski}
#'    \item{Xie_Beni}\item{GDI11}\item{GDI12}\item{GDI13}\item{GDI21}\item{GDI22}
#'    \item{GDI23}\item{GDI31}\item{GDI32}\item{GDI33}\item{GDI41}\item{GDI42}
#'    \item{GDI43}\item{GDI51}\item{GDI52}\item{GDI53}
#' }
#' @export
#'
SamplingLouvain <- function(dataMatrix = dataMatrix,
                            outlier_kmeans = NULL,
                            knn_range = c(3:70),
                            iter = 30,
                            compute_index = "Calinski_Harabasz",
                            cores = 1,
                            seed = NULL){

  int_index = compute_index
  knn_n = knn_range
  if(!is.null(seed)){
    set.seed(seed)
    seeds = sample(1:1e8, iter)
  }else{
    seeds = sample(1:1e8, iter)
  }
  ds_kmeans_df_keep = outlier_kmeans$ds_kmeans_df_keep
  ds_kmeans = outlier_kmeans$ds_kmeans

  # prepare iterations function for parallel
  Iterations = function(s){

    cat("\n\n*** iter: ",s,"***\n")
    set.seed(seeds[s])
    sampledf = strata(ds_kmeans_df_keep,
                      stratanames = "cluster",
                      size = rep(1,nrow(ds_kmeans$centers)),
                      method = "srswor",
                      description = FALSE)
    sampledf$cellname = rownames(ds_kmeans_df_keep)[sampledf$ID_unit]
    sample_centers = dataMatrix[sampledf$cellname,]
    rownames(sample_centers) = 1:nrow(sample_centers)
    sample_centers = apply(sample_centers,2,function(x){
      as.numeric(as.character(x))
    })

    cat("\n - knn + louvain \n")
    d_mick = parDist(x = sample_centers, method = "euclidean", diag = TRUE, upper = TRUE, threads = 1)
    louvain_list = lapply(knn_n, function(k){
      if( k < nrow(as.matrix(d_mick))){
        knn = dbscan::kNN(d_mick, k = k)
        # knn = dbscan::sNN(d_mick, k = k)
        knn_edge = reshape2::melt(knn$id)
        knn_graph = graph_from_edgelist(as.matrix(knn_edge[,c(1,3)]), directed = FALSE)
        c_louvain = cluster_louvain(knn_graph)
        c_louvain
      }else{
        NULL
      }
    })
    names(louvain_list) = paste0("k", knn_n)

    cat("\n - calculate internal indicators \n")
    cluster_meta_small = do.call(cbind, lapply(louvain_list, function(l){ l$membership }) )
    rownames(sample_centers) = paste0("cell",1:nrow(sample_centers))
    rownames(cluster_meta_small) = paste0("cell",1:nrow(sample_centers))

    best_df = FindBestCluster(data = sample_centers,
                              meta_data = cluster_meta_small,
                              methods = int_index,
                              cores = 1)
    result = list(knn_n = knn_n,
                  sampledf = sampledf,
                  louvain_list = louvain_list,
                  int_result = best_df)
    return(result)
  }

  cl = makeCluster(cores, outfile = "Sampling_Log.txt")
  clusterExport(cl = cl,
                varlist = c("FindBestCluster", "strata", "parDist", "graph_from_edgelist", "cluster_louvain",
                            "detectCores", "intCriteria", "bestCriterion"))
  clusterExport(cl = cl, varlist = c("ds_kmeans_df_keep","ds_kmeans","seed","dataMatrix","knn_n","int_index"), envir=environment())
  # clusterEvalQ(cl, library(clusterCrit))
  cat("\n/// sampling \n")
  sampling_list = parLapply(cl, 1:iter, Iterations)
  stopCluster(cl)

  int_best = IntIndexCounting(sampling_list)
  result = list(sampling_list = sampling_list,
                int_best = int_best)

  return(result)
}

#' Count internal index results
#'
#' @description
#' It is a built-in function of SamplingLouvain for voting statistics of internal indicator results
#'
#' @param sampling_list result from SamplingLouvain function
#'
#' @export
#'
IntIndexCounting <- function(sampling_list = NULL){

  cat("\n/// obtain the best number of neighbors \n")
  int_best = do.call(cbind,lapply(sampling_list, function(x){
    as.character(x$int_result$best)
  }) ) %>% as.data.frame
  int_best$best = apply(int_best, 1, function(y){names(which.max(table(y)))})
  rownames(int_best) = rownames(sampling_list[[1]]$int_result)

  return(int_best)
}

#' Find best cluster result from internal indicators
#'
#' @Description:
#' This function is a built-in function of SamplingLouvain,
#' which used to find the best clustering result by
#' counting the results of internal indicators
#'
#' @param data row is cell, col is feature
#' @param meta_data columns represent different clustering results
#' @param methods internal index
#' @param cores  the number of threads
#'
#' @export
#'
FindBestCluster <- function(data = NULL,
                            meta_data = NULL,
                            methods = c("all"),
                            cores = 1){

  max_cores = detectCores()-2
  if(is.null(data)){
    stop("Data is null.")
  }
  if(is.null(meta_data)){
    stop("Metadata is null.")
  }
  if(is.null(rownames(meta_data))){
    stop("meta_data missing row name.")
  }
  # if(!all.equal(rownames(meta_data),rownames(data))){
  #   stop("The cell name of the data is inconsistent with the row name (cell name) of the metadata.")
  # }


  orig_internal_index = c("Ball_Hall", "Banfeld_Raftery", "C_index", "Calinski_Harabasz", "Davies_Bouldin",
                          "Det_Ratio", "Dunn", "Gamma", "G_plus", "GDI11", "GDI12", "GDI13", "GDI21", "GDI22",
                          "GDI23", "GDI31", "GDI32", "GDI33", "GDI41", "GDI42", "GDI43", "GDI51", "GDI52",
                          "GDI53", "Ksq_DetW", "Log_Det_Ratio", "Log_SS_Ratio", "McClain_Rao", "PBM", "Point_Biserial",
                          "Ray_Turi", "Ratkowsky_Lance", "Scott_Symons", "SD_Scat", "SD_Dis", "S_Dbw", "Silhouette",
                          "Tau", "Trace_W", "Trace_WiB", "Wemmert_Gancarski", "Xie_Beni")

  tmp = setdiff(methods, orig_internal_index)
  if(length(tmp)==0){
    crit = methods
  }else if((tmp=="all")){
    crit = orig_internal_index
  }else{
    stop(paste0("There is no such indicator: ", tmp))
  }

  # get cluster_name
  if(is.null(colnames(meta_data))){
    cluster_name = paste0("cluster_result_",1:ncol(meta_data))
  } else{
    cluster_name = colnames(meta_data)
  }

  # compute parallel
  if(cores == 1){
    result = do.call(cbind, lapply(as.list(cluster_name), function(i){
      unlist(intCriteria(data, as.integer(as.numeric(as.character(meta_data[,i]))), crit = crit))
    }))
  }else if(cores < max_cores){
    cl = makeCluster(cores)
    registerDoParallel(cl)
    result = foreach(i = cluster_name, .combine = 'cbind') %dopar% {
      unlist(clusterCrit::intCriteria(data, as.integer(as.numeric(as.character(meta_data[,i]))), crit = crit)) }
    stopCluster(cl)
  }else{
    cat("exceed max cores, use cores = 1")
    result = do.call(cbind, lapply(as.list(cluster_name), function(i){
      unlist(intCriteria(data, as.integer(as.numeric(as.character(meta_data[,i]))), crit = crit))
    }))
  }
  colnames(result) = cluster_name

  # find best cluster number
  bestCluster <- mapply(function(x, y){
    val = result[x,][!(is.na(result[x,])|is.infinite(result[x,]))]
    index = bestCriterion(val, y)
    names(val)[index]},
    as.list(1:length(crit)), crit)

  # analysis result
  result = data.frame(result,
                      methods = rownames(result),
                      best = bestCluster,
                      stringsAsFactors = FALSE)
  rownames(result) = crit

  return(result)

}
