# HCMMCNVs

HCMMCNVs is a browser-based software for detecting copy number variants (CNVs) using whole exome sequencing (WES) technology in combining multiple processed bam inputs with different disease or tumor types.

This repository is a fork of the `lunching/HCMM_CNVs` repository.  This repository contains code that allows HCMM-CNVs to be run as a command line tool, rather than a shiny app.  
Please cite the original authors if you use HCMM-CNVs in your work.  

Song, C., Su, S.C., Huo, Z., Vural, S., Galvin, J.E. and Chang, L.C., 2021. 
HCMMCNVs: hierarchical clustering mixture model of copy number variants detection using whole exome sequencing technology. 
Bioinformatics. https://doi.org/10.1093/bioinformatics/btab183

#### Requirement
* R >= 4.0.1
* `Rsamtools` package

## Running the HCMM-CNVs wrapper

Make sure R is loaded in your environment and then run the following in the folder in which you would like your results to be stored.  
```
bedfile="/path/to/HCMM_CNVs_wrapper/Toy_example/Demo.bed"
chr="19"
bam_directory="/path/to/HCMM_CNVs_wrapper/Toy_example/"
min_coverage="10"
filename="Test"
n_cluster="3"
code_dir="/path/to/HCMM_CNVs_wrapper/"
ploidy_file="1"

Rscript --vanilla /path/to/HCMM_CNVs_wrapper/HCMM_CNVs_wrapper.R $code_dir \
$bedfile $chr $bam_directory $min_coverage $filename $n_cluster $ploidy_file
```
