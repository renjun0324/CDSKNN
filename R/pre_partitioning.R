
#' K-means clustering based on outlier detection
#'
#' @description
#' Get the result of K-means area division after outlier detection.
#'
#' @param dat expression matrix, row is cell, col is feature
#' @param partition_count the number of partitions
#' @param num_init the proportion of outliers that need to be removed
#' @param batch_size the number of threads
#' @param max_iters the maximum iterations for kmeans
#' @param outlier_det whether running outlier detection
#' @param outlier_methods outlier detection methods, md or ed.
#' "md" utilizes the Mahalanobis distance for detection,
#' while "ed" treats points that are far from the center proportionally as outliers.
#' @param outlier_q the proportion of outliers you want to filter.
#' When outlier_methods is ed, this parameter takes effect.
#' @param min_n the minimum number of data points for outlier detection.
#' @param seed random seed
#'
#' @return A list that include K-means result, division result, and center matrix was removed outliers
#' @importFrom pbapply pblapply
#' @export
#'
pre_partitioning <- function(dat = NULL,
                             partition_count = NULL,
                             batch_size = min(500, nrow(dat)),
                             num_init = 10,
                             max_iters = 1000,
                             outlier_det = TRUE,
                             outlier_methods = "md",
                             outlier_q = 0.2,
                             min_n = 200,
                             seed = 723){

  k = partition_count

  #######################  1. kmeans
  del_dat <- t(dat)
  set.seed(seed)
  mb_r <- mbkmeans(del_dat,
                   clusters = k,
                   batch_size = as.integer(batch_size),
                   num_init = num_init,
                   max_iters = max_iters )
  dat_l = tapply(1:nrow(dat),
                 mb_r$Clusters,
                 function(x,y){ dat[x,,drop = FALSE] })
  cluster <- sapply(1:length(dat_l), function(i){
    dat = dat_l[[i]]
    tmp = rep(i, nrow(dat))
    names(tmp) = rownames(dat)
    return(tmp)
  }) %>% unlist
  centers <- do.call(rbind, lapply(dat_l, function(dat){
    apply(dat, 2, mean)
  }) )
  kmeans_result <- list(cluster = cluster,
                        centers = centers,
                        data_list = dat_l)

  #######################  2. outlier detect
  if(outlier_det){
    outlier_s = pblapply( 1:length(dat_l), function(i){
      dat = dat_l[[i]]
      outlier = tryCatch({
        if(outlier_methods == "md"){ OutlierDetMD(dat, min_n = min_n) }
        if(outlier_methods == "ed"){ OutlierDetED(dat, q = outlier_q, min_n = min_n) }
      }, error = function(e){
        cat("error in detect outlier", i, "\n")
      })
      keep = setdiff(rownames(dat), outlier)
      center = apply(dat[keep,,drop=F], 2, mean)
      return(list(outlier = outlier,
                  keep = keep,
                  center = center))
    })
    outlier_point = sapply(outlier_s, function(x) x$outlier)
    keep_point = sapply(outlier_s, function(x) x$keep)
    keep_centers = do.call(rbind, lapply(outlier_s, function(x) x$center))
    rownames(keep_centers) = 1:nrow(keep_centers)
    outlier_result <- list(outlier = outlier_point,
                           stay = keep_point,
                           centers = keep_centers,
                           cluster = kmeans_result$cluster[unlist(keep_point)])
  }else{
    outlier_result <- NULL
  }

  #######################  3. combine result
  if(outlier_det){
    cluster_df <- data.frame(row.names = names(outlier_result$cluster),
                             cellname = names(outlier_result$cluster),
                             cluster = outlier_result$cluster )
  }else{
    cluster_df <- data.frame(row.names = names(kmeans_result$cluster),
                             cellname = names(kmeans_result$cluster),
                             cluster = kmeans_result$cluster )
  }
  result = list(kmeans_result = kmeans_result,
                outlier_result = outlier_result ,
                cluster_df = cluster_df)

  #######################  4. command
  command = GetCommand()
  command$params = command$params[-1] %>% unlist
  result$command = command
  return(result)
}

#' Outlier Detection
#'
#' @description
#' This function is used to detect outliers of any matrix.
#'
#' @param m matrix, row is cells
#' @param q the proportion of outliers you want to filter
#' @param min_n the minimum number of data points for outlier detection.
#' @param cores the number of threads
#'
#' @return Character vector of outlier label
#' @export
#'
OutlierDetED <- function(m,
                         q = 0.2,
                         min_n = 5){

  # calculate the sum of the squares from the center to each point
  if(is.null(rownames(m))){
    names = 1:nrow(m)
  }else{
    names = rownames(m)
  }
  max_cores = detectCores()-2
  if(!is.matrix(m)){
    m = as.matrix(m)
  }

  if(q==0){
    final = NULL
  }else{
    if(nrow(m) <= min_n){
      final = NULL
    }else{
      m_m = base::colMeans(m)
      m2 = m - t(m_m %*% matrix(1, 1, nrow(m)))
      m3 = rowMeans(m2^2)
      max = quantile(m3, 1-q)
      ind = which(m3 > max)
      final = names[ind]
    }
  }
  return(final)
}


#' Outlier Detection
#'
#' @description
#' This function is used to detect outliers of any matrix,
#' which comes from mvoutlier,
#' please check [github] (https://github.com/cran/mvoutlier) for more function details.
#'
#' @param m matrix, row is cells
#' @param min_n the minimum number of cells required
#' @param alpha parameter of covMcd
#' @param pvalue significance level.
#'
#' @importFrom robustbase covMcd
#'
#' @return Character vector of outlier label
#' @export
#'
OutlierDetMD <- function(m,
                         min_n = 5,
                         alpha = 1,
                         pvalue = 0.05){

  if(nrow(m) <= min_n){
    final = NULL
  }

  if(nrow(m)<=ncol(m)){
    final = NULL
  }else if(nrow(m)==ncol(m)+1){
    final = NULL
  }else{

    if(is.null(rownames(m))){
      names = 1:nrow(m)
    }else{
      names = rownames(m)
    }
    covr = covMcd(m, alpha=alpha)
    dist = mahalanobis(m, center=covr$center, cov=covr$cov)
    p_value = 1 - pchisq(dist, ncol(m)-1)
    outl = which(p_value < pvalue)

    final = names[outl]
    if(length(final)==0){
      final = NULL
    }
  }

  return(final)
}
