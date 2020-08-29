require('rjson')
require('tidyverse')



grab_sample_df <- function(my_sample, study_acc){
  sample_acc <- paste(study_acc, my_sample$assay$name, sep="_")
  #print(sample_acc)
  t_df <- transpose(my_sample$characteristic)
  t_df$category[sapply(t_df$category, is.null)] <- ""
  t_df$value[sapply(t_df$value, is.null)] <- ""
  n <- length(my_sample$characteristic)
  df <- tibble(
    "study_acc"=rep(study_acc, n),
    "sample_acc"=rep(sample_acc, n),
    "key"=unlist(t_df$category),
    "value"=unlist(t_df$value))
  return(df)
}

grab_study_df <- function(study_acc){
  print(study_acc)
  exp_dat <- fromJSON(file=sprintf("data/ae_out/%s.json", study_acc) )
  samples <- exp_dat$experiment$sample
  df <- samples %>% map(~grab_sample_df(.x, study_acc)) %>% bind_rows()
  return(df)
}

ae_studies <- read_tsv("data/array_express_studies.tsv", col_names=FALSE) 
ae_attrib <- ae_studies %>% pull(X1) %>% map(grab_study_df) %>% bind_rows()
ae_attrib %>% write_csv("data/ae_attributes.csv")