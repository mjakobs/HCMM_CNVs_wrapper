library(Rsamtools)

HCMM_CNVs<- function(Cov_matrix, n_cluster, bed_file_sorted, sample_names, filename, progress){
  # Hierarchical Clustering
  h_all<- c()
  n_cluster_all<- c()
  for(i in 1: ncol(Cov_matrix)){
    progress$inc((1/ncol(Cov_matrix))*0.4 , detail = paste("Running Hierarchical Clustering for region ", i, sep=""))
    hc_tmp<- hclust(dist(Cov_matrix[,i]/median(Cov_matrix[,i])))
    h_tmp<- 0
    n_cluster_tmp <- length(unique(cutree(hc_tmp, h = h_tmp)))
    while(n_cluster_tmp>3){
      h_tmp<- h_tmp + 0.05
      n_cluster_tmp <- length(unique(cutree(hc_tmp, h = h_tmp)))
    }
    h_all[i]<- h_tmp
    n_cluster_all[i]<- n_cluster_tmp
  }
  
  mean_height<- mean(h_all)
  trim_cov_mean_hc_ref<- c()
  
  for(i in 1:ncol(Cov_matrix)){
    hc_tmp<- hclust(dist(Cov_matrix[,i]/median(Cov_matrix[,i])))
    cluster_tmp<- cutree(hc_tmp, h = mean_height) 
    cluster_max_tmp<- which.max(table(cluster_tmp))
    trim_cov_mean_hc_ref[i]<- mean(Cov_matrix[which(cluster_tmp == cluster_max_tmp),i])
  }
  matrix_adj_Cov_hc_log2_ratio<- log(t(Cov_matrix)/trim_cov_mean_hc_ref, base = 2)
  
  ## Run Circular Binary Segmentation
  maploc_center_tmp<- (bed_file_sorted[,2] + bed_file_sorted[,3])/2
  CBS_all<- list()
  
  for(i in 1: ncol(matrix_adj_Cov_hc_log2_ratio)){
    progress$inc((1/ncol(matrix_adj_Cov_hc_log2_ratio))*0.4 , detail = paste("Running CBS for sample ", i, sep=""))
    log_ratio_adj_cov_hc_tmp<- matrix_adj_Cov_hc_log2_ratio[, i]
    CNA_data_hc_tmp<- CNA(log_ratio_adj_cov_hc_tmp, chrom = bed_file_sorted[,1], maploc = maploc_center_tmp, data.type = "logratio", sampleid = sample_names[i])
    smoothed_CNA_hc_object <- smooth.CNA(CNA_data_hc_tmp)
    segment_smoothed_CNA_hc_object <- segment(smoothed_CNA_hc_object, verbose=1)
    CBS_all[[i]]<- segment_smoothed_CNA_hc_object
  }

  progress$inc(0.1, detail = paste("Saving result into RData"))
  names(CBS_all) <- sample_names
  save(CBS_all, file = paste("CBS_", filename, ".RData", sep=""))
  test<- paste("There are ", ncol(Cov_matrix), " target regions.", sep="")
  return(list(test = test))
}




