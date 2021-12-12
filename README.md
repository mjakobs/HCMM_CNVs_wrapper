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

You will need to provide the following information:
* `bedfile`: A bed file 
* `chr`: The chromosome that you would like to investigate, e.g. `"19"`.  Alternatively, if `chr` is set to `"all"`, all chromosomes present in the bed file will be analysed.  
* `bam_directory`: A `.txt` file containing the full path to the bam files that you would like to investigate.  One bam file per row - see the `input_bams.txt` file in the `Toy_example` folder for proper formatting
* `min_coverage`: Minimum mean coverage.  The default value is `"10"`.  HCMM-CNVs will filter regions below the threhold in all samples. 
* `filename`: The name of your analysis
* `n_clusters`: The number of clusters for the hierarchical clustering step.  This should be a value from 2-4.  The default value is `"3"`
* `code_dir`: The full path to where you have installed `HCMM_CNVs_wrapper`
* `ploidy_file`: This should be `"1"` if you do not have ploidy information, or the path to a file containing your ploidy information. 

Make sure R is loaded in your environment and then run the following in the folder in which you would like your results to be stored.  
```
bedfile="/path/to/HCMM_CNVs_wrapper/Toy_example/Demo.bed"
chr="19"
bam_directory="/path/to/HCMM_CNVs_wrapper/Toy_example/input_bams.txt"
min_coverage="10"
filename="Test"
n_cluster="3"
code_dir="/path/to/HCMM_CNVs_wrapper/"
ploidy_file="1"

Rscript --vanilla /path/to/HCMM_CNVs_wrapper/HCMM_CNVs_wrapper.R $code_dir \
$bedfile $chr $bam_directory $min_coverage $filename $n_cluster $ploidy_file
```
