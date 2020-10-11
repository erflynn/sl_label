require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]


my_files <- list.files(sprintf("data/rnaseq/%s/01_study_mat/", prefix))
res <- lapply(my_files, function(my_f){
       dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s", prefix, my_f));
       return(list("samples"=colnames(dat), "nrow"=nrow(dat)))})
save(res, file=sprintf("data/03_qc/%s_counts_rna.RData", prefix))