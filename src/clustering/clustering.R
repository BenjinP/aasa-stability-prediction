

#### Library loading ####
source("C:/Users/Benjamin/Documents/Repos/proyecto_titulacion/ProteinStabilityPrediction/R/library_loading.R")

datasets_names = c('p1STN', 'p4LYZ', 'p1BPI', 'HLYZ')

dataset_name = 'HLYZ'
dataset_plot_name = '4LYZ'
data_url = 'https://raw.githubusercontent.com/Naio/aasa-stability-prediction/master/data/processed/original/'
dataset_url = paste(data_url, dataset_name,'.csv', sep='')
protein_data = read.csv(file=dataset_url)

#row.names(protein_data) = paste(case_when(
#  protein_data$stability < 0 ~ "(-)",
#  protein_data$stability == 0 ~ "(=)", 
#  protein_data$stability > 0 ~ "(+)"),
#  protein_data$id,
#  sep=""
#)


#Changes row names to show the ID and the stability of the mutant, so in the clusters visualization
#that information appears
row.names(protein_data) <- paste("(",protein_data$stability,") ", protein_data$id, sep="")
protein_data$id = NULL


protein_data <- protein_data[order(protein_data$stability),]
descriptors <- protein_data[,-c(1,2)]
standarized_descriptors = as.data.frame(scale(descriptors))




#### Data preprocessing ####
library(dplyr)



#For discretizacion, different breaks and labels are define
labeled_protein_data <- mutate(protein_data, class = case_when(
  stability < 0 ~ "Negative variation",
  stability == 0 ~ "No variation", 
  stability > 0 ~ "Positive variation"))

labeled_protein_data$class <- as.factor(labeled_protein_data$class)




#### Clustering ####
#Libraries
library(cluster)
library(gridExtra)
library(factoextra)
library(amap)
#Settings for every algorithm
clustering_seed <- 100
clustering_data <- standarized_descriptors
distance_method <- "euclidean"
dist_mat <- Dist(clustering_data, method = distance_method)
min_n_clusters <- 2
max_n_clusters <- 12
RESULTS_DIR = "results/"


#The matrix indicates how many instances of each type
#(Positive Variation, Negative Variation, No variation) are in each cluster. 
create_conf_mat <- function(labeled_data ,clusters) {
  #Making pseudo confusion matrix
  conf_mat <- table(labeled_data[,"class"], clusters)
  conf_mat
}

#Needs the algorithm nane
save_conf_mat <- function(file_path, conf_mat) {
  #Opens a png file to write on it
  png(file_path,height = 100, width = 330)
  #Outputs the table to the png file
  grid.table(conf_mat)
  #Closes the png file
  dev.off()
}

##### KMeans ####

#This is a wrapper function for the Kmeans function, which allows to have
#control of the distance function the algorithm is using.
kmeans_wrapper <- function(data, k, distance_method, nstart, iter.max) {
  set.seed(clustering_seed)
  Kmeans(x=data, centers = k, method = distance_method, 
         nstart = nstart, iter.max = iter.max)
}

#This function calculates and plot the silhouette coefficient 
#for clustering configurations producer with k values between 1 and 10.
#Silhouette coefficient evaluates how "good" are the clusters for a 
#given k value.
kmeans_silhouette_plot <- fviz_nbclust(clustering_data, kmeans_wrapper, 
                                  method = "silhouette", 
                                  distance_method = distance_method, 
                                  nstart = 1000, iter.max = 500)

for (k in 2:12) {
  set.seed(100)
  Kmeans.out <- Kmeans(x=clustering_data, centers = k, method = distance_method, nstart = 1000, iter.max = 500)
  cluster_plot <- fviz_cluster(Kmeans.out, data = standarized_descriptors, repel=TRUE) + ggtitle(paste("K-Means clustering with K = ",k," for Dataset ", dataset_plot_name, sep=""))
  ggsave(plot = cluster_plot, paste(RESULTS_DIR, "Kmeans_",k,".png",sep=""), width = 7, height = 6)
  
  #Making pseudo confusion matrix
  conf_mat <- create_conf_mat(labeled_protein_data, Kmeans.out$cluster)
  #Saving confusion matrix
  save_conf_mat(file_path = paste(RESULTS_DIR, "Kmeans_",k,"_conf_",dataset_name,".png",sep=""),
                conf_mat = conf_mat)
}



#### Hierarchical ####

#The fviz_nbclust will pass the data, not the distance metric, so
#as a workaround, we pass it as an aditional argument
hcut_wrapper <- function(data, k, linkage_method, distance_mat) {
  #Deterministic, so no seed needed.
  hcut(distance_mat, k, hc_func = "hclust", hc_method = linkage_method)
}

#These are the linkage methods used. average is UPGMA, and mcquitty is WPGMA
linkage_methods <- c("complete")

for (linkage_method in linkage_methods) {
  hclust_result <- hclust(dist_mat, method=linkage_method)
  
  results_folder <- paste(RESULTS_DIR, "Hclust", sep="")
  #Creates a folder to save the clustering results
  dir.create(path = results_folder)
  
  #Cutting the dendrogram.
  for(k in min_n_clusters:max_n_clusters) {
    #Plotting the dendrograms
    hclust_plot <- fviz_dend(hclust_result, k=k, cex = 0.7, repel=T, lwd = 0.8,
                             rect = T, rect_fill = T,
                             main = paste("Hierarchical clustering with k =", k,
                                          "for Dataset", dataset_plot_name ))
    
    ggsave(plot = hclust_plot, paste(results_folder, "/HClust_",dataset_name,"_k",k,".png",sep=""), 
                                      width = 16, height = 6)
    
    hcut_results <- hcut_wrapper(clustering_data, k, linkage_method, dist_mat)
    #Making pseudo confusion matrix
    hcut_results$cluster
    conf_mat <- create_conf_mat(labeled_protein_data, hcut_results$cluster)
    #Saving confusion matrix
    
    save_conf_mat(file_path = paste(results_folder, "/CONF_",
                                    dataset_name,"_HClust_k",k,".png",sep=""),
                  conf_mat = conf_mat)
    
  }
  
}


#### Visualizing distances matrices ####
library(heatmaply)


#The path of the folder should be specified
#A suffix for every plot name.
calculate_dist_matrices <- function(data, path, plot_suffix) {
   
  #Visualizing distances
  methods <-c("euclidean", "maximum", "manhattan", "pearson",
             "abspearson", "correlation", "abscorrelation", "spearman", "kendall")
 
  names(methods) <- c("Euclidean", "Maximum", "Manhattan", 
                     "Not Centered Pearson", "Absolute Pearson", "Centered Pearson", 
                     "Absolute correlation", "Spearman", "Kendall")
   
  for(method.name in names(methods)) {
    
    #Calculates distance matrix between observations in datasets
    distance_matrix <- as.matrix(Dist(data, method=methods[method.name]))  
    
    #Plotting and saving to directory
    heatmaply(distance_matrix,
              scale="none", dendrogram = "none",
              main=paste("Distance matrix",  method.name, "method", "for Dataset", dataset_plot_name),
              revC = T, column_text_angle = 90,
              cexCol = 0.6, cexRow = 0.6,
              #With these parameters, heatmaply function can save to disk
              file = paste(path, methods[method.name], plot_suffix,".png", sep=""),
              width = 1000, height = 1000
              )
    
  }
}

library(heatmaply)
library(orca)# To export ggplots
calculate_dist_matrices(standarized_descriptors,"results/", dataset_name)





