library(Rsamtools)

data_preproc<- function(bed_file, chr, bam_dir, min_cov, filename){
  # bam file directory
  #progress$inc(0.1, detail = paste("Checking bam files"))
  print(paste("Checking bam files"))
  text_bam_dir<- paste(bam_dir)
  # Check how many bam files in bam file directory
  # bam_dir_tmp <-  dir(bam_dir)[grep(".bam",dir(bam_dir), fixed = T)]
  # bam_dir_tmp2 <- bam_dir_tmp[grep(".bam",sapply(1: length(bam_dir_tmp), function(x) substr(bam_dir_tmp[x], nchar(bam_dir_tmp[x])-3, nchar(bam_dir_tmp[x])) ))]
  bam_dir_tmp <- read.table(bam_dir, header = FALSE)
  bam_dir_tmp2 <- bam_dir_tmp[,1]
  text_bam_numbers<- paste("There are ", length(bam_dir_tmp2), " bam files.", sep="" )

  # Read in bed file
  bed_file <- read.table(bed_file, sep = "\t", header = T)
  # Check given bed files
  # bed_file <- bed_file[bed_file[,1]==chr, ]
  bed_file_sorted <- bed_file[order(bed_file[,1],bed_file[,2]), ]
  text_bed_regions<- paste("There are ", nrow(bed_file_sorted), " regions in exome ", chr, ".", sep="")
  # Generate ID
  bam_dir_tmp3 <- sub(".*/", "", bam_dir_tmp2)
  id_tmp <- gsub(".bam", "", bam_dir_tmp3)
  # Generate the range from bed file
  ref_tmp <- IRanges(start = bed_file_sorted[, 2], end = bed_file_sorted[, 3])
  # Calculate the coverage of all targed regions for all individuals and total mapped read (TMR)
  start_tmp <- start(ref_tmp)[1]
  end_tmp <- end(ref_tmp)[length(ref_tmp)]
  which <- IRangesList(quack = IRanges(start_tmp - 10000, end_tmp + 10000))
  names(which) <- as.character(chr)
  what <- c("pos", "mapq", "qwidth")
  flag <- scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE, isFirstMateRead = TRUE)
  param <- ScanBamParam(which = which, what = what, flag = flag)
  sample_names<- sub(".bam", "", bam_dir_tmp3)
  Cov_matrix <- matrix(NA, nrow = length(ref_tmp), ncol = length(bam_dir_tmp2))
  TMR <- rep(0, length(bam_dir_tmp2))
  
  for(i in 1: length(bam_dir_tmp2)){
   # progress$inc((1/length(bam_dir_tmp2))*0.6 , detail = paste("Calculating coverage and total mapped read of sample ", i, sep=""))
    print(paste("Calculating coverage and total mapped read of sample ", i, sep=""))
    bam <- scanBam(paste(bam_dir_tmp2[i]), param = param)[[1]]
    #readlength[i] <- round(mean(bam[["qwidth"]]))
    irang <- IRanges(bam[["pos"]], width = bam[["qwidth"]])
    Cov_matrix[,i] <- countOverlaps(ref_tmp, irang)
    TMR[i] <- sum(idxstatsBam(paste(bam_dir_tmp2[i]))[,3])
  }
  
  ##----- Data pre-processing -----##
  # 1. Remove duplicated regions
  #progress$inc(0.05, detail = paste("Removing duplicated regions..."))
  print(paste("Removing duplicated regions..."))
  index_duplicated_region<- intersect(which(duplicated(bed_file_sorted[,2])), which(duplicated(bed_file_sorted[,3])))
  
  if(length(index_duplicated_region)>0){
    Cov_matrix<- Cov_matrix[-index_duplicated_region, ]
    #GC_all<- GC_all[-index_duplicated_region]
    bed_file_sorted<- bed_file_sorted[-index_duplicated_region, ]
  }
  # 2. Calculate the mean coverages for all region accross all cell lines 
  #progress$inc(0.05, detail = paste("Calculating the mean coverages for all regions..."))
  print(paste("Calculating the mean coverages for all regions..."))
  mean_cov<- apply(Cov_matrix, 1, mean) # mean coverages for all regions
  index_mean_cov_10<- which(mean_cov < 10) # Index for the region with mean coverage from all regions across all cell lines less than 10
  
  if(length(index_mean_cov_10)>0){
    #GC_all<- GC_all[-index_mean_cov_10]
    Cov_matrix<- Cov_matrix[-index_mean_cov_10,]
    bed_file_sorted<- bed_file_sorted[-index_mean_cov_10, ]
  }
  # 3. Calculate the adjusted coverage 
  #progress$inc(0.1, detail = paste("Calculating adjusted coverage for all regions..."))
  print(paste("Calculating adjusted coverage for all regions..."))
  mean_LS<- mean(TMR) # Mean of the library size
  lib_ratio<- mean_LS/TMR # ratio of the library size of each samples  
  matrix_adj_Cov_rm_duplicated<-  t(Cov_matrix + 1)*lib_ratio # +1 coverage to aviod NA 
  
  text_cov_summary<- paste("After data processing, there are ", ncol(matrix_adj_Cov_rm_duplicated), " target regions.", sep="")
  #progress$inc(0.1, detail = paste("Generating coverage matrix..."))
  print(paste("Generating coverage matrix..."))
  save(sample_names, bed_file_sorted, matrix_adj_Cov_rm_duplicated, file = paste("Cov_matrix_", filename, ".RData", sep=""))
  text_cov_rdata<- paste("Coverage results were saved to ",getwd(),"Cov_matrix_", filename, ".RData.",  sep="")
  
  #text1<- head(bed_file)
  #text2<- paste("chr", chr, sep="")
  #text3<- length(dir(bam_dir))
  #text4<- paste("Minimum coverage: ", min_cov, sep="")
  #return(list(text1 = text1, text2 = text2, text3 = text3, text4 = text4))
  return(list(text_bam_dir = text_bam_dir, text_bam_numbers = text_bam_numbers, text_bed_regions = text_bed_regions, text_cov_summary = text_cov_summary, text_cov_rdata = text_cov_rdata))
}
