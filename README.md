# CDSKNN Clustering 

Large-scale single-cell clustering algorithm based on K-means and optimal KNN graph structure

## Installation

```py
pip install faiss-cpu --user
```

```r
install.packages("devtools")
devtools::install_github("renjun0324/CDSKNN")
```

## Quick Start

```r
data(pca_result)
data(cellinfo)

result = cdsknn(dat = pca_result, 
                partition_count = 200, 
                batch_size = 500,
                num_init = 10,
                max_iters = 10,
                outlier_det = TRUE,
                outlier_methods = "md",
                outlier_q = 0.2, 
                min_n = 200,
                cluster_method = "louvain",
                resolution = 1,
                knn_range = c(3:20), 
                iter = 50,
                is_weight = TRUE,
                assess_index = "Calinski_Harabasz", 
                res_range = seq(0.2, 2, 0.1),
                new_cluster_method = "louvain",
                seed = 723)


# ARI result
library(aricode)
ARI(result$cluster_result$cluster[rownames(cellinfo)],cellinfo$celltype)
```

## Steps

```r
# outlier detect
outlier_kmeans = pre_partitioning(dat = pca_result,
                                  partition_count = 300,
                                  batch_size = 300,
                                  outlier_det = T,
                                  outlier_methods = "md",
                                  outlier_q = 0.2,
                                  min_n = 200,
                                  num_init = 10,
                                  max_iters = 10,
                                  seed = 723)
      
# random sampling and choose optmial KKN graph structure                      
sampling_result = resampling(dat = pca_result, 
                             outlier_kmeans = outlier_kmeans, 
                             resolution = 1, 
                             cluster_method = "louvain", 
                             knn_range = c(3:20), 
                             iter = 20, 
                             is_weight = TRUE, 
                             seed = 723)
                                  
# get final clustering result                                 
new_r = new_clustering(sampling_result = sampling_result,
                       outlier_kmeans = outlier_kmeans,
                       assess_index = "Calinski_Harabasz",
                       cluster_method = "louvain",
                       res_range = seq(0.2, 2, 0.1),
                       is_weight = TRUE,
                       seed = 723)
                       
result <- list(outlier_kmeans = outlier_kmeans,
               sampling_result = sampling_result,
               cluster_result = new_r)

# ARI result
library(aricode)
ARI(result$cluster_result$cluster[rownames(cellinfo)],cellinfo$celltype)
```
