# KKL Clustering

Large-scale single-cell clustering algorithm based on K-means and optimal KNN graph structure

## Installation

```r
install.packages("devtools")
devtools::install_github("renjun0324/KKLClustering")
```

## Quick Start

```r
data(pca_result)
data(cellinfo)

result = kkl(pca_result,
             outlier_q = 0.1,
             down_n = 300,
             knn_range = 5:70,
             iter = 100,
             compute_index =  c("Davies_Bouldin","Calinski_Harabasz"),
             assess_index = "Davies_Bouldin",
             cores = 1,
             seed = 723)

# ARI result
library(aricode)
ARI(result$cluster_df$cluster,cellinfo$celltype)
```

## Steps

```r
# outlier detect
outlier_kmeans = OutlierKmeans(dataMatrix = pca_result,
                               outlier_q = 0.1,
                               down_n = 300,
                               cores = 1,
                               seed = 723)
      
# random sampling and choose optmial KKN graph structure                      
sampling_result = SamplingLouvain(dataMatrix = pca_result,
                                  outlier_kmeans = outlier_kmeans,
                                  knn_range = 5:70,
                                  iter = 100,
                                  compute_index = c("Davies_Bouldin","Calinski_Harabasz"),
                                  cores = 10,
                                  seed = 723)
                                  
# get final clustering result                                 
new_louvain = NewLouvain(sampling_result = sampling_result,
                         outlier_kmeans = outlier_kmeans,
                         assess_index = "Davies_Bouldin",
                         cores = 1)
                       
result = data.frame(row.names = rownames(new_louvain$cluster_meta),
                    name = rownames(new_louvain$cluster_meta),
                    celltype = cellinfo$celltype,
                    cluster = as.factor(new_louvain$cluster_meta$cluster),
                    stringsAsFactors = FALSE)

# ARI result
library(aricode)
ARI(result$cluster, result$celltype)
# 0.8486607                               
```
