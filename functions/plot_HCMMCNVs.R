plot_HCMMCNVs<- function(result, sample_id){
  CBS_tmp<- result[[match(sample_id, names(result))]]
  return(list(CBS_tmp = CBS_tmp))
}