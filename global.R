## Check packages ##
if (!require('DT')) install.packages("DT")
if (!require('shiny')) install.packages("shiny") 
if (!require('shinythemes')) install.packages("shinythemes") 
if (!require('shinyFiles')) install.packages("shinyFiles") 

if (!require('DNAcopy')){
  source("https://bioconductor.org/biocLite.R")
  biocLite("DNAcopy")
}
 
if (!require('Rsamtools')){
  source("https://bioconductor.org/biocLite.R")
  biocLite("Rsamtools")
}

if (!require('httpuv')){
  install.packages("httpuv")
}

## R functions ##
source("functions/data_preproc.R")
source("functions/HCMM_CNVs.R")
source("functions/plot_HCMMCNVs.R")
source("functions/Mixture_Model.R")
