# Rscript to run HCMM-CNVs from the command line
source("functions/data_preproc.R")
source("functions/HCMM_CNVs.R")
source("functions/Mixture_Model.R")
source("functions/plot_HCMMCNVs.R")

# expected input preproc: bed_file, chr, bam_dir, min_cov, filename -- removed progress
# expected output preproc: sample_names, bed_file_sorted, matrix_adj_Cov_rm_duplicated
# expected input HCMMCNVs: Cov_matrix (matrix_adj_Cov_rm_duplicated), n_cluster, bed_file_sorted (bed_file_sorted), sample_names (sample_names), filename, progress, ploidy = "1", ploidy_path = ""
# expected output HCMMCNVs: CBS_all
# expected input Mixture_Model: means, markers, tau, n.states = 3, common.sigma2 = FALSE, tol = 1e-6

# expected input plot_HCMMCNVs:result, sample_id

vars.tmp <- commandArgs()
vars <- vars.tmp[length(vars.tmp)]
split.vars <- unlist(strsplit(vars," "))
###  Setting input paths for normalized read count and experimental design ###
bed_file <- split.vars[1]
chr <- split.vars[2]
bam_dir <- split.vars[3]
min_cov <- split.vars[4]
filename <- split.vars[5]
n_cluster <- split.vars[6]

######################
# Run pre-processing #
######################
data_preproc(bed_file = bed_file, chr = chr, bam_dir = bam_dir, min_cov = min_cov, filename = filename, )
