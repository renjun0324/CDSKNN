
#' K-means clustering based on outlier detection
#'
#' @description
#' Get the result of K-means area division after outlier detection.
#'
#' @param dataMatrix row is cell, col is feature
#' @param down_n downsampling number
#' @param outlier_q the proportion of outliers that need to be removed
#' @param cores the number of threads
#' @param seed random seed
#'
#' @return A list that include K-means result, division result, and center matrix was removed outliers
#' @export
#'
OutlierKmeans <- function(dataMatrix = NULL,
                          down_n = 300,
                          outlier_q = 0.2,
                          cores = 1,
                          seed = 723){

  cat("\n/// Kmeans with", down_n, "centers \n")
  set.seed(seed)
  ds_kmeans = kmeans(dataMatrix, centers = down_n, iter.max = 1000)
  if(is.null(names(ds_kmeans$cluster))){
    names(ds_kmeans$cluster) = rownames(dataMatrix)
  }
  ds_kmeans_df = data.frame(row.names = names(ds_kmeans$cluster),
                            cellname =  names(ds_kmeans$cluster),
                            cluster = ds_kmeans$cluster,
                            stringsAsFactors = FALSE)
  ds_kmeans_df = ds_kmeans_df[order(ds_kmeans_df$cluster),]

  cat("\n/// outlier detect \n")
  # debug
  # r=tapply(1:nrow(dataMatrix),
  #          ds_kmeans$cluster,
  #          function(x,y){
  #            x
  #        })
  ds_cluster_list = tapply(1:nrow(dataMatrix),
                           ds_kmeans$cluster,
                           function(x,y){
                             # cat(y,"\n")
                             if(length(x)==1){
                               # cat(y,"\n")
                               tmp = as.matrix(dataMatrix[x,])
                               # tmp = as.matrix(dataMatrix[x,])
                               rownames(tmp) = rownames(dataMatrix)[x]
                               tmp
                             }else{ dataMatrix[x,] }})

  cl = makeCluster(cores, outfile = "OutlierKmeans_Log.txt")
  clusterExport(cl = cl, varlist = c("OutlierDet", "outlier_q","ds_cluster_list"), envir=environment())
  #clusterExport(cl = cl, varlist = c("ds_cluster_list"), envir=environment())
  clusterEvalQ(cl, library(doParallel))
  outlier_list_2 = parLapply(cl,
                             1:length(ds_cluster_list),
                             function(j){
                               outlier = tryCatch({OutlierDet(ds_cluster_list[[j]], q = outlier_q)},
                                                  error = function(e){
                                                    cat("error in detect outlier", j, "\n") })
                               keep = setdiff(rownames(ds_cluster_list[[j]]), outlier)
                               center = tryCatch({
                                 if(length(keep)==1){
                                   tt = data.frame(ds_cluster_list[[j]][keep,])
                                   colnames(tt) = keep
                                   tt = t(tt)
                                 }else{
                                   tt = data.frame(ds_cluster_list[[j]][keep,])
                                 }
                                 apply(tt, 2, mean)},
                                 error = function(e){
                                   cat("error in center without outlier", j, "\n")
                                 })
                               return(list(outlier = outlier,
                                           keep = keep,
                                           center = center))

                             })
  stopCluster(cl)

  #cat("\n/// debug \n")
  # get result
  outlier_list = lapply(outlier_list_2, function(x) x$outlier)
  keep_list = lapply(outlier_list_2, function(x) x$keep)
  ds_kmeans_df_keep = ds_kmeans_df[unlist(keep_list),]
  ds_kmeans_centers = do.call(rbind, lapply(outlier_list_2, function(x) x$center))
  rownames(ds_kmeans_centers) = 1:nrow(ds_kmeans_centers)

  outlier_kmeans = list(ds_kmeans = ds_kmeans,
                        ds_cluster_list = ds_cluster_list,
                        outlier = outlier_list,
                        ds_kmeans_df = ds_kmeans_df,
                        ds_kmeans_df_keep = ds_kmeans_df_keep,
                        ds_kmeans_new_centers = ds_kmeans_centers)

  return(outlier_kmeans)
}



#' Outlier Detection
#'
#' @description
#' This function is used to detect outliers of any matrix.
#'
#' @param m matrix, row is cells
#' @param q the proportion of outliers you want to filter
#' @param min_n the minimum number of cells required
#' @param cores the number of threads
#'
#' @return Character vector of outlier label
#' @export
#'
OutlierDet <- function(m,
                       q = 0.1,
                       min_n = 3,
                       cores = 1){

  # for each point, calculate the sum of the squares from
  # the original center to the center without the point
  max_cores = detectCores()-2
  if(!is.matrix(m)){
    m = as.matrix(m)
  }

  if(nrow(m) <= min_n){
    final = NULL
  }else{
    m_m = colMeans(m)

    if(cores == 1){
      m2 = foreach(x = 1:nrow(m), .combine = 'rbind') %do% { colMeans(m[-x,]) - m_m }
    }else if(cores < max_cores){
      cl = makeCluster(cores)
      registerDoParallel(cl)
      m2 = foreach(x = 1:nrow(m), .combine = 'rbind') %dopar% { colMeans(m[-x,]) - m_m }
      stopCluster(cl)
    }else{
      cat("exceed max cores, use cores = 1")
      m2 = foreach(x = 1:nrow(m), .combine = 'rbind') %do% { colMeans(m[-x,]) - m_m }
    }

    m3 = rowMeans(m2^2)
    max = quantile(m3, 1-q)
    ind = which(m3 > max)

    final = rownames(m)[ind]
  }

  return(final)
}
