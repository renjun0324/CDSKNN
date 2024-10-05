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

result = cdsknn(dat = pca, 
                partition_count = 300, 
                batch_size = 200,
                num_init = 10,
                max_iters = 100,
                outlier_det = TRUE,
                outlier_methods = "md",
                outlier_q = 0.2, 
                min_n = 200,
                cluster_method = "louvain",
                resolution = 1,
                knn_range = c(3:10), 
                iter = 20,
                is_weight = TRUE,
                assess_index = "Calinski_Harabasz", 
                res_range = seq(0.2, 2, 0.2),
                new_cluster_method = "louvain",
                python_path = "/usr/bin/python3",
                seed = 723)


# ARI result
library(aricode)
ARI(result$cluster_result$cluster[rownames(water)],water$true_cluster)
```

## Steps

```r
# outlier detect
outlier_kmeans = pre_partitioning(dat = pca,
                                  partition_count = 300,
                                  batch_size = 200,
                                  outlier_det = T,
                                  outlier_methods = "md",
                                  outlier_q = 0.2,
                                  min_n = 200,
                                  num_init = 10,
                                  max_iters = 100,
                                  seed = 723)
      
# random sampling and choose optmial KKN graph structure                      
sampling_result = resampling(dat = pca, 
                             outlier_kmeans = outlier_kmeans, 
                             resolution = 1, 
                             cluster_method = "louvain", 
                             knn_range = c(3:10), 
                             iter = 20, 
                             is_weight = TRUE, 
                             python_path = "/usr/bin/python3",
                             seed = 723)
                                  
# get final clustering result                                 
new_r = new_clustering(sampling_result = sampling_result,
                       outlier_kmeans = outlier_kmeans,
                       assess_index = "Calinski_Harabasz",
                       cluster_method = "louvain",
                       res_range = seq(0.2, 2, 0.2),
                       is_weight = TRUE,
                       # python_path = "/usr/bin/python3",
                       seed = 723)
                       
result <- list(outlier_kmeans = outlier_kmeans,
               sampling_result = sampling_result,
               cluster_result = new_r)

# ARI result
library(aricode)
ARI(result$cluster_result$cluster[rownames(water)],water$true_cluster)
```

## Plots

```r

library(ggplot2)
library(ggsci)
library(patchwork)

umap = umap::umap(pca)
df = data.frame(umap$layout,
                cluster = result$cluster_result$cluster[rownames(water)],
                true_cluster = water$true_cluster)
df$cluster = as.factor(df$cluster)

p1 <- ggplot(df, aes(x = X1, y = X2, color = cluster)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_d3("category20") +
  labs(x = "UMAP-1", y = "UMAP-2", title = "CDSKNN Clustering", color = "") +
  theme(panel.background = element_rect(fill='transparent', color="black"),
        legend.key = element_rect(fill = "white"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.key.size = unit(0.4, "cm") )

p2 <- ggplot(df, aes(x = X1, y = X2, color = true_cluster)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_d3("category20") +
  labs(x = "UMAP-1", y = "UMAP-2", title = "True Cluster", color = "") +
  theme(panel.background = element_rect(fill='transparent', color="black"),
        legend.key = element_rect(fill = "white"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        legend.key.size = unit(0.4, "cm") )
        
p = p1 + p2
ggsave("umap.png", p, width = 8, height = 3)
```
