

#' Random sampling to do KNN-Louvain/leiden clustering
#'
#' @description
#' This function is used to randomly sample each cluster multiple times
#' to select the most powerful KNN graph structure.
#'
#' @param dat expression matrix with cell * feature
#' @param outlier_kmeans the result from pre_partitioning function
#' @param knn_range the range of the number of neighbors in the KNN graph structure
#' @param cluster_method louvain or leiden in resampling
#' @param resolution clustering parameters settings in resampling
#' @param iter the number of iterations in resampling
#' @param is_weight Whether to use distance weights when constructing the KNN graph
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
resampling <- function(dat = dat,
                       outlier_kmeans = NULL,
                       knn_range = c(3:70),
                       cluster_method = "louvain",
                       resolution = 1,
                       iter = 30,
                       is_weight = TRUE,
                       python_path = "/usr/bin/python3",
                       seed = 723){

  cluster_df = outlier_kmeans$cluster_df

  #######################  1. set random seeds
  down_number = nrow(outlier_kmeans$kmeans_result$centers)
  if(!is.null(seed)){
    set.seed(seed)
    seeds = sample(1:1e8, iter)
  }else{
    seeds = sample(1:1e8, iter)
  }

  #######################  2. ampling + create graph
  Iterations = function(s, resolution){

    # 1). -- sampling
    strata = rep(1, down_number)
    names(strata) = 1:down_number
    set.seed(seeds[s])
    sample_v = stratified(cluster_df, "cluster", strata)
    sample_centers = dat[sample_v$cellname,] %>% as.matrix
    # sample_centers = cellnames[sample_v$cellname] %>% as.matrix

    # 2). -- knn
    g_list = lapply(knn_range, function(k){
      # cat(k, "\n")
      g = CreateKNN(sample_centers, k, is_weight = TRUE, python_path = python_path)
      return(g)
    })
    names(g_list) = paste0("k",knn_range)

    # 3). -- cluster
    cl_list = lapply(knn_range, function(k){
      # cat(k,"\n")
      g = g_list[[paste0("k",k)]]
      if(cluster_method == "louvain"){
        set.seed(seed)
        cl = cluster_louvain(g, weights = E(g)$weight, resolution = resolution)$membership
      }
      if(cluster_method == "leiden"){
        set.seed(seed)
        cl = cluster_leiden(g, weights = E(g)$weight, resolution_parameter = resolution)$membership
      }
      return(cl)
    })
    cl_df = do.call(cbind, cl_list)
    colnames(cl_df) = paste0("k", knn_range)
    return(list(g_list = g_list, cluster = cl_df))
  }
  cluster_iters = pblapply(1:iter, function(i){ Iterations(i, resolution = resolution)  })

  #######################  3. RMI
  cat(knn_range, "\n")
  k_rmi_list = pblapply(knn_range, function(k){
    t = do.call(cbind, lapply(1:iter, function(itr){
      cluster_iters[[itr]]$cluster[,paste0("k",k)]
    }))
    m = matrix(NA, nr = iter, nc = iter)
    for(i in 1:(iter-1)){
      for(j in (i+1):iter){
        # m[i,j] = reduced_mutual_information(t[,i],t[,j], method = "approximation1", normalized = TRUE)
        m[i,j] = clustAnalytics::reduced_mutual_information(t[,i],t[,j], method = "approximation2", normalized = TRUE)
      }
    }
    upper_t = m[upper.tri(m)]
    return(upper_t)
  })

  #######################  4. output best_k
  k_rmi = lapply(k_rmi_list, function(m){
    median(m)
  }) %>% unlist
  best_k = knn_range[which.max(k_rmi)]
  cat("best k: ", best_k, "\n")

  #######################  5. get sampling result
  sampling_result <- list(sampling_cluster = cluster_iters,
                          neighbors_RMI = k_rmi,
                          best_neighbors = best_k)

  #######################  6. command
  command = GetCommand()
  command$params = command$params[-c(1:2)]
  sampling_result$command = command

  return(sampling_result)
}

#' faiss_ann
#'
#' @description
#' The invocation interface of faiss in Python.
#'
#' @param dat used to create knn graph
#' @param k the result from mbkmeans function
#' @param python_path python path to be used
#' @importFrom reticulate import
#' @export
#'
faiss_ann <- function(dat, k, python_path){

  if(!is.null(python_path)){
    use_python(python_path)
  }

  k = k + 1
  faiss = import('faiss')
  np = import("numpy", convert=FALSE)

  x = np_array(dat, dtype = "float32")
  index = faiss$IndexFlatL2(as.integer(ncol(dat)))
  index$add(x)
  y = index$search(x, as.integer(k))
  y[[1]] = y[[1]][,-1]
  y[[2]] = y[[2]][,-1] + 1
  return(y)
}


#' Create knn-graph
#'
#' @description
#' This function is used to running Louvain clustering
#'
#' @param dat row is cell, col is feature
#' @param k the number of nearest neighbors
#' @param is_weight whether to use distance weights when constructing the KNN graph
#' @param python_path python path to be used
#' @export
#'
CreateKNN <- function(dat,
                      k,
                      is_weight = FALSE,
                      python_path = "/usr/bin/python3"){

  faiss_r <- faiss_ann(dat, k = k+1, python_path)

  if(is_weight){
    edge = reshape2::melt(faiss_r[[2]])
    weight = reshape2::melt(faiss_r[[1]])
    g = graph_from_edgelist(as.matrix(edge[,c(1,3)]), directed = FALSE)
    E(g)$weight = weight[,3]
  }else{
    edge = reshape2::melt(faiss_r[[2]])
    g = graph_from_edgelist(as.matrix(edge[,c(1,3)]), directed = FALSE)
  }
  return(g)
}

#' stratified sampling
#'
#' @description
#' stratified sampling function from github：
#' https://gist.github.com/mrdwab/6424112
#'
#' @param df The input data.frame
#' @param group A character vector of the column or columns that make up the "strata".
#' @param size The desired sample size.
#‘ If size is a value less than 1, a proportionate sample is taken from each stratum.
#' If size is a single integer of 1 or more, that number of samples is taken from each stratum.
#' If size is a vector of integers, the specified number of samples is taken for each stratum. It is recommended that you use a named vector. For example, if you have two strata, "A" and "B", and you wanted 5 samples from "A" and 10 from "B", you would enter size = c(A = 5, B = 10).
#' @param select This allows you to subset the groups in the sampling process. This is a list. For instance, if your group variable was "Group", and it contained three strata, "A", "B", and "C", but you only wanted to sample from "A" and "C", you can use select = list(Group = c("A", "C")).
#' @param replace For sampling with replacement.
#' @export
#'
stratified <- function(df,
                       group,
                       size,
                       select = NULL,
                       replace = FALSE,
                       bothSets = FALSE) {
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)

  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}

