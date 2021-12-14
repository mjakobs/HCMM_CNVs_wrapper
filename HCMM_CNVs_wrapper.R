# Rscript to run HCMM-CNVs from the command line
library(Rsamtools)
library(DNAcopy)
options(bitmapType='cairo')

# expected input preproc: bed_file, chr, bam_dir, min_cov, filename -- removed progress
# expected output preproc: sample_names, bed_file_sorted, matrix_adj_Cov_rm_duplicated
# expected input HCMMCNVs: Cov_matrix (matrix_adj_Cov_rm_duplicated), n_cluster, bed_file_sorted (bed_file_sorted), sample_names (sample_names), filename, progress, ploidy = "1", ploidy_path = ""
# expected output HCMMCNVs: CBS_all
# expected input Mixture_Model: means, markers, tau, n.states = 3, common.sigma2 = FALSE, tol = 1e-6

# expected input plot_HCMMCNVs:result, sample_id

# getting command line arguments
vars <- commandArgs(trailingOnly=TRUE)
###  Setting input variables from jobscript ###
code_dir <- vars[1]
bed_file <- vars[2]
chr <- vars[3]
bam_dir <- vars[4]
min_cov <- vars[5]
filename <- vars[6]
n_cluster <- vars[7]
ploidy_path <- vars[8]

#set working directory to code directory so you can source the functions
source(paste(code_dir,"functions/data_preproc.R", sep = ""))
source(paste(code_dir,"functions/HCMM_CNVs.R", sep = ""))
source(paste(code_dir,"functions/Mixture_Model.R", sep = ""))
source(paste(code_dir,"functions/plot_HCMMCNVs.R", sep = ""))

if(chr != "all"){
  print("Running single chromosome analysis")
  chr_tmp <- chr
  ######################
  # Run pre-processing #
  ######################
  data_preproc(bed_file = bed_file, chr = chr_tmp, bam_dir = bam_dir, min_cov = min_cov, filename = filename)
  
  print("Finished pre-processing!")
  
  ######################
  #    Run HCMM-CNVs   #
  ######################
  
  ploidy_option <- NA
  if(paste(ploidy_path) == "1"){
    ploidy_option <- "ploidy = \"1\""
  }else{
    ploidy_option <- "ploidy_path = ploidy_path"
  }
  paste(ploidy_option)
  load(file = paste("Cov_matrix_chr_",chr_tmp,"_", filename, ".RData", sep=""))
  
  HCMM_CNVs(Cov_matrix = matrix_adj_Cov_rm_duplicated, n_cluster = n_cluster, 
            bed_file_sorted = bed_file_sorted, sample_names = sample_names, 
            filename = filename, paste(ploidy_option), chr = chr_tmp)
  
  ######################
  # Run plot HCMM-CNVs #
  ######################
  
  # don't actually need to run the R function provided, use the DNAcopy plot instead
  
  load(file = paste("CBS_chr_",chr,"_", filename, ".RData", sep=""))
  
  for(i in 1:length(names(CBS_all))){
    png(paste("CNVs_chr_",chr,"_", names(CBS_all)[i], ".png", sep=""))
    plot(CBS_all[[i]])
    dev.off()
  }
  
  # save text copy number output
  for(t in 1:length(names(CBS_all))){
    sample_tmp <- names(CBS_all)[t]
    write.table(CBS_all[[t]][[2]], paste(sample_tmp,"_chr_",chr,"_copy_number.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
}else if(chr == "all"){
  print("Analysing all chromosomes")
  
  bed <- read.table(bed_file, sep = "\t", header = T)
  chr_list <- levels(as.factor(bed[,1]))
  chr_list <- chr_list[nchar(chr_list)<=5]
  chr_list <- chr_list[-grep("M",chr_list)]
  
  # make output table for each sample
  bam_names <- read.table(bam_dir, header = F, sep = "\t")
  bam_names[,1] <- gsub(".*/","",bam_names[,1])
  bam_names[,1] <- gsub(".bam","",bam_names[,1])
  
  
  cn_data <- as.data.frame(matrix(nrow = 1, ncol = 6))
  colnames(cn_data) <- c("ID", "chrom", "loc.start",  "loc.end", "num.mark", "seg.mean")
  
  list_df <- lapply(bam_names[,1], function(x) cbind(bam_names[x,1],cn_data))
  names(list_df) <- paste(bam_names[,1],"_CN", sep = "")
  list_df <- lapply(names(list_df),function(x) list_df[[x]][,-1])
  names(list_df) <- paste(bam_names[,1],"_CN", sep = "")
  list_df <- lapply(names(list_df),function(x) list_df[[x]][-1,])
  names(list_df) <- paste(bam_names[,1],"_CN", sep = "")
  
  for(c in 1:length(chr_list)){
    chr_tmp <- chr_list[c]
    ######################
    # Run pre-processing #
    ######################
    data_preproc(bed_file = bed_file, chr = chr_tmp, bam_dir = bam_dir, min_cov = min_cov, filename = filename)
    
    print("Finished pre-processing!")
    
    ######################
    #    Run HCMM-CNVs   #
    ######################
    
    ploidy_option <- NA
    if(paste(ploidy_path) == "1"){
      ploidy_option <- "ploidy = \"1\""
    }else{
      ploidy_option <- "ploidy_path = ploidy_path"
    }
    paste(ploidy_option)
    load(file = paste("Cov_matrix_chr_", chr_tmp,"_", filename, ".RData", sep=""))
    
    HCMM_CNVs(Cov_matrix = matrix_adj_Cov_rm_duplicated, n_cluster = n_cluster, 
              bed_file_sorted = bed_file_sorted, sample_names = sample_names, 
              filename = filename, paste(ploidy_option), chr = chr_tmp)
    
    ######################
    # Run plot HCMM-CNVs #
    ######################
    
    # don't actually need to run the R function provided, use the DNAcopy plot instead
    
    load(file = paste("CBS_chr_",chr_tmp,"_", filename, ".RData", sep=""))
    
    # save copy number figure
    for(i in 1:length(names(CBS_all))){
      png(paste("CNVs_chr_",chr_tmp,"_", names(CBS_all)[i], ".png", sep=""))
      plot(CBS_all[[i]])
      dev.off()
    }
    # save copy number text
    for(s in 1:length(names(list_df))){
      sample_tmp <- names(list_df)[s]
      sample_tmp <- gsub("_CN","",sample_tmp)
      list_df[[s]] <- rbind(list_df[[s]],CBS_all[[sample_tmp]][[2]])
      
    }
  }
  for(t in 1:length(names(list_df))){
    sample_tmp <- names(list_df)[t]
    sample_tmp <- gsub("_CN","",sample_tmp)
      write.table(list_df[[t]], paste(sample_tmp,"_copy_number.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# Clean up intermediate files
sapply(list.files(pattern = ".RData"), unlink)

