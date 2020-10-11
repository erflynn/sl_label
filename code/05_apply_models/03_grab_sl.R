# grab the predicted 
require('tidyverse')
prefix <- "human"
data_type <- "microarray"
ds <- "sex"

label_dir <- "data/08_model_out"
my_f <- list.files(path=label_dir, pattern=sprintf("%s_%s_%s*", prefix, data_type, ds))
stopifnot(length(my_f) >0) # // TODO - check size depending on dataset

sl_df <- do.call(rbind, lapply(my_f, function(f){
  read_csv(sprintf("%s/%s", label_dir, f))
}))

write_csv(sl_df, sprintf("data/09_model_summary/%s_%s_%s_labels.csv", prefix, data_type, ds))