# figure out which RNA-seq files are present

require('tidyverse')

options(stringsAsFactors=FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

load(sprintf("data/03_qc/%s_counts_rna.RData", prefix))

summary(sapply(res, function(x) x$nrow))

samples <- unique(unlist(sapply(res, function(x) x$samples)))
samples2 <- samples[samples!="gene_name"]
sample_metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix) )
sample_metadata2 <- sample_metadata %>% mutate(present=(acc %in% samples2)) 
table(sample_metadata2$present) # 122864 out of 229789
# mouse:  125652 out of 359142
# human: 6252 out of 615+6252
mapping <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample.csv", prefix))
mapping2 <- mapping %>% mutate(present=(sample_acc %in% samples2))
mapping2 %>% write_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample2.csv", prefix))
sample_metadata2 %>% write.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata2.csv", prefix))

mapping2 %>% filter(present) %>% select(study_acc) %>% unique() %>% count() # 6290 human, 6499 mouse
mapping2 %>% select(study_acc) %>% unique() %>% count() # 6457 human, 6651 mouse, 396 rat -- all present

# check #2
list.studies <- unique(mapping$study_acc)
num_samples <- lapply(list.studies, function(study){
  list.samples <- list.files(path=sprintf("data/rnaseq/%s/00_infiles/%s/",  prefix, study),
                             pattern="*.sf")
  return(length(list.samples))
})

num_per_study <- mapping %>% group_by(study_acc) %>% count()

df <- data.frame(cbind("study"=unlist(list.studies), "num_s"=unlist(num_samples)))
num_files_present <- sum(unlist(num_samples)) # mouse: 128079, human: 122885

counts_f <- num_per_study %>% 
  left_join(df, by=c("study_acc"="study")) %>% 
  rename(num_samples=n, num_files=num_s ) %>%
  mutate(num_files=as.numeric(as.character(num_files))) %>%
  mutate(num_missing=num_samples-num_files)

counts_f %>% arrange(desc(num_missing))

counts_f %>% filter(num_files==100 & num_missing > 100) %>% 
  write_csv(sprintf("%s_over_100_missing.csv", prefix))


