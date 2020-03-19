# figure out which RNA-seq files are present

load(sprintf("data/%s_counts_rna.RData", prefix))
summary(sapply(res, function(x) x$nrow))

samples <- unique(unlist(sapply(res, function(x) x$samples)))
samples2 <- samples[samples!="gene_name"]
sample_metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix) )
sample_metadata2 <- sample_metadata %>% mutate(present=(acc %in% samples2)) 
table(sample_metadata2$present) # 122864 out of 229789
# mouse:  125652 out of 359142
mapping <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample.csv", prefix))
mapping2 <- mapping %>% mutate(present=(sample_acc %in% samples2))
mapping2 %>% write_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample2.csv", prefix))
sample_metadata2 %>% write.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata2.csv", prefix))

mapping2 %>% filter(present) %>% select(study_acc) %>% unique() %>% count() # 6290 human, 6499 mouse
mapping2 %>% select(study_acc) %>% unique() %>% count() # 6457 human, 6651 mouse

