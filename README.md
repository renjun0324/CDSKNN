# KKL Clustering

Large-scale single-cell clustering algorithm based on K-means and optimal KNN graph structure

## Installation

```r
install.packages("devtools")
devtools::install_github("renjun0324/KKLClustering")
```

## Usage

```r
data(pca_result)
data(cellinfo)

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
                                  int_index = c("Davies_Bouldin","Calinski_Harabasz"),
                                  cores = 10,
                                  seed = 723)
                                  
# get final clustering result                                 
new_louvain = NewLouvain(sampling_result = sampling_result,
                         outlier_kmeans = outlier_kmeans,
                         index = "Davies_Bouldin",
                         cores = 1)
result = data.frame(row.names = rownames(new_louvain$cluster_meta),
                    name = rownames(new_louvain$cluster_meta),
                    cluster = as.factor(new_louvain$cluster_meta$cluster),
                    stringsAsFactors = FALSE)
                               
```
