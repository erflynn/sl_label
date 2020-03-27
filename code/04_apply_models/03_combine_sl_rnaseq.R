# put together the RNAseq output labels

require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

my_path <- sprintf("data/rnaseq/%s/04_model_out/", prefix)
list_f <- list.files( path=my_path, pattern=sprintf("%s_rnaseq_sl*", prefix))

f_df <- do.call(rbind, lapply(list_f, function(f) read_csv(sprintf("%s/%s", my_path, f))))
f_df %>% write_csv(sprintf("data/02_imputed_labels/%s_rnaseq_sl.csv", prefix))