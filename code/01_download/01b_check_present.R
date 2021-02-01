library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

setwd(sprintf("data/03_expression/rnaseq/%s", prefix))
rp <- list.files(pattern="^[E|D|S]RP*")
res <- lapply(rp, function(i) list.files(i, pattern="*.sf"))
list_present <- str_replace_all(unlist(res), "_quant.sf", "")
data.frame(rp[which(lapply(res, length)==0)]) %>% 
  write_csv(sprintf("../%s_studies_missing_samples.csv",prefix))
data.frame(list_present) %>% 
  write_csv(sprintf("../%s_samples_present.csv", prefix))