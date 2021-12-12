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
  
}else if(chr == "all"){
  print("Analysing all chromosomes")
  
  bed <- read.table(bed_file, sep = "\t", header = T)
  chr_list <- levels(as.factor(bed[,1]))
  for(i in 1:length(chr_list)){
    chr_tmp <- chr_list[i]
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
    
    for(i in 1:length(names(CBS_all))){
      png(paste("CNVs_chr_",chr_tmp,"_", names(CBS_all)[i], ".png", sep=""))
      plot(CBS_all[[i]])
      dev.off()
    }
    
  }
  
}

