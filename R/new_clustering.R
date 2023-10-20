
#' Find best cluster result from internal indicators
#'
#' @Description:
#' This function is a built-in function of pre_partitioning,
#' which used to find the best clustering result by
#' counting the results of internal indicators
#'
#' @param data row is cell, col is feature
#' @param meta_data columns represent different clustering results
#' @param methods internal index
#'
#' @export
#'
FindBestCluster <- function(data = NULL,
                            meta_data = NULL,
                            methods = c("all")){

  if(is.null(data)){
    stop("Data is null.")
  }
  if(is.null(meta_data)){
    stop("Metadata is null.")
  }
  if(is.null(rownames(meta_data))){
    stop("meta_data missing row name.")
  }

  orig_internal_index = c("Ball_Hall", "Banfeld_Raftery", "C_index", "Calinski_Harabasz", "Davies_Bouldin",
                          "Det_Ratio", "Dunn",  "Gamma", "G_plus", "Ksq_DetW", "Log_SS_Ratio", "McClain_Rao",
                          "PBM", "Point_Biserial", "Ratkowsky_Lance", "Ray_Turi",  "Scott_Symons",
                          "SD_Scat", "SD_Dis", "S_Dbw", "Silhouette", "Tau",
                          "Trace_W", "Trace_WiB", "Wemmert_Gancarski", "Xie_Beni"
                          # "GDI11", "GDI12", "GDI13", "GDI21", "GDI22", "GDI53",
                          # "GDI23", "GDI31", "GDI32", "GDI33", "GDI41", "GDI42", "GDI43", "GDI51", "GDI52",
                          #  "Log_Det_Ratio",
  )

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

  # compute
  result = do.call(cbind, lapply(as.list(cluster_name), function(i){
    unlist(intCriteria(data, as.integer(as.numeric(as.character(meta_data[,i]))),
                       crit = crit))
  }))
  colnames(result) = cluster_name

  # find best cluster number
  bestCluster <- mapply(function(x, y){
    val = result[x,][!(is.na(result[x,])|is.infinite(result[x,]))]
    # val = val[which(val!=0)]
    # val = val[which(val!=1)]
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

#' Louvain/leiden clustering with optimal KNN structure
#'
#' @description
#'
#' This function is used to run Louvain clustering with optimal KNN graph structure
#'
#' @param sampling_result sampling result from pre_partitioning function
#' @param outlier_kmeans outlier detection result from OutlierKmeans function
#' @param assess_index evaluation index used to select the optimal KNN graph structure
#' @param cluster_method louvain or leiden
#' @param res_range resolution range in clustering
#' @param seed random seed
#'
#' @export
#'
new_clustering <- function(sampling_result = NULL,
                           outlier_kmeans = NULL,
                           assess_index = "Calinski_Harabasz",
                           cluster_method = "louvain",
                           res_range = seq(0.2,3,0.2),
                           is_weight = TRUE,
                           seed = 723){

  kmeans_result = outlier_kmeans$kmeans_result

  #######################  0. get best k
  best_k = sampling_result$best_neighbors
  iter = length(sampling_result$sampling_cluster)

  #######################  1. clustering with different resolution
  pblapply(res_range, function(res){
    # cat(res, "\n")
    sapply(1:iter, function(i){
      # cat(i, "\n")
      g = sampling_result$sampling_cluster[[i]]$g_list[[paste0("k",best_k)]]
      if(cluster_method == "louvain"){
        set.seed(seed)
        cl = cluster_louvain(g, weights = E(g)$weight, resolution = res)$membership
      }
      if(cluster_method == "leiden"){
        set.seed(seed)
        cl = cluster_leiden(g, weights = E(g)$weight, resolution_parameter = res)$membership
      }
      return(cl)
    }) -> meta
    dimnames(meta) = list(rownames(outlier_kmeans$kmeans_result$centers),1:iter)
    return(meta)
  }) -> meta_list
  names(meta_list) = paste0("res", res_range)

  #######################  2. ch result
  pblapply(meta_list, function(meta){
    best_df = FindBestCluster(data = kmeans_result$centers,
                              meta_data = meta,
                              methods = c("Calinski_Harabasz") )
    return(best_df)
  }) -> ch_list
  ch = do.call(rbind, ch_list)[,1:iter]
  ch_means = rowMeans(ch, na.rm = F)
  best_res = res_range[which.max(ch_means)]

  #######################  3. best cluster
  g = CreateKNN(kmeans_result$centers, best_k, is_weight = is_weight)
  if(cluster_method == "louvain"){
    set.seed(seed)
    best_cluster = cluster_louvain(g, weights = E(g)$weight, resolution = best_res)$membership
  }
  if(cluster_method == "leiden"){
    set.seed(seed)
    best_cluster = cluster_leiden(g, weights = E(g)$weight, resolution_parameter = best_res)$membership
  }
  names(best_cluster) = 1:length(best_cluster)
  final = best_cluster[kmeans_result$cluster]
  names(final) = names(kmeans_result$cluster)
  cat("best res:",best_res, "\n")

  #######################  4. cluster_result
  cluster_result = list(res_result = ch,
                        best_res = best_res,
                        cluster = final[outlier_kmeans$cluster_df$cellname])

  #######################  6. command
  command = GetCommand()
  command$params = command$params[-c(1:2)]
  cluster_result$command = command

  return(cluster_result)
}
