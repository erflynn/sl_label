# Try alternate RNA-seq models
library(tidyverse)
load("data/07_model_dat/human_rnaseq_sex_train_mat.RData") # X_train, X_test, Y_train, Y_test

list_rows <- tibble("sample_acc"=c(rownames(X_train), rownames(X_test)))
sample_meta <- read_csv("data/01_sample_lists/list_samples_w_qc.csv")
samples_tt <- sample_meta %>% semi_join(list_rows)
samp_to_study <- read_csv("data/01_sample_lists/sample_to_study.csv")

# 0. Remove scRNA-seq
# TODO update this
sample_tt_v2 <- samples_tt %>% filter(!scrna, present, read_count > 1000000)
list_train_test <- samp_to_study %>% semi_join(sample_tt_v2, by="sample_acc")
list_train_test %>% distinct(study_acc, sample_acc) %>% write_csv("data/04_test_sets/human_rnaseq_tt_v2.csv")


sample_to_remove <- sample_to_study %>% semi_join(to_remove, by="study_acc") %>% distinct(sample_acc)
train_to_remove <- intersect(rownames(X_train), sample_to_remove$sample_acc)
train_idx <- which(rownames(X_train) %in% sample_to_remove$sample_acc)
test_to_remove <- intersect(rownames(X_test), sample_to_remove$sample_acc)
test_idx <- which(rownames(X_test) %in% sample_to_remove$sample_acc)
X_train2 <- X_train[-train_idx,]
Y_train2 <- Y_train[-train_idx]
X_test2 <- X_test[-test_idx,]
Y_test2 <- Y_test[-test_idx]



# 1. filter for different counts reads
samp_meta <- read_csv("data/sample_metadata_filt.csv", col_types="cccccdldcc")
samp_reads <- samp_meta %>% filter(organism=="human", data_type=="rnaseq") %>%
  select(sample_acc, num_reads)
train_reads <- samp_reads %>% filter(sample_acc %in% rownames(X_train2))
test_reads <- samp_reads %>% filter(sample_acc %in% rownames(X_test2))

# 1mil
CUTOFF <- 1000000
train_exc <- train_reads %>% filter(num_reads<CUTOFF) %>% pull(sample_acc)
test_exc <- test_reads %>% filter(num_reads<CUTOFF) %>% pull(sample_acc)
train_idx2 <- which(rownames(X_train2) %in% train_exc)
test_idx2 <- which(rownames(X_test2) %in% test_exc)
X_train3 <- X_train2[-train_idx2,]
Y_train3 <- Y_train2[-train_idx2]
X_test3 <- X_test2[-test_idx2,]
Y_test3 <- Y_test2[-test_idx2]



# 2. check on PAR region mapping
load("data/00_reference/genes/xy_genes_annot_par_biomart.RData") # --> human_rnaseq
my_transcripts <- colnames(X_train3) # 2906
hr2 <- human_rnaseq %>% filter(ensembl_transcript_id %in% my_transcripts) # 2887
hr2 %>% group_by(chromosome_name, region) %>% count()

# do female samples have any counts mapped to Y chromsome??
f_samp <- X_train3[Y_train3==0,]
m_samp <- X_train3[Y_train3==1,]

f_samp_y <- f_samp[,(hr2 %>% filter(chromosome_name=="Y") %>% pull(ensembl_transcript_id))]
head(f_samp_y)

m_samp_y <- m_samp[,(hr2 %>% filter(chromosome_name=="Y") %>% pull(ensembl_transcript_id))]
head(m_samp_y)

# TODO: need to look at actual read counts
human_escape <- read_csv("data/00_reference/genes/aggreg_human_escape.csv")

# REMOVE PAR
par_remove <- which(colnames(X_train) %in% (hr2 %>% filter(region!="non-PAR")%>% pull(ensembl_transcript_id)))
X_train4 <- X_train3[,-par_remove]
X_test4 <- X_test3[,-par_remove]
data.frame(X_train4) %>% write_csv("data/rnaseq_xtrain_filt.csv")
data.frame(Y_train3) %>% write_csv("data/rnaseq_ytrain_filt.csv")

x_to_study <- sample_to_study %>% 
  filter(sample_acc %in% rownames(X_train4)) %>%
  group_by(sample_acc) %>% 
  summarize(study_acc=paste(unique(study_acc), collapse=";")) # none get collapsed
df <- data.frame(x_to_study)
rownames(df) <- df$sample_acc
df2 <- df[rownames(X_train4),]


df2 %>% write_csv("data/rnaseq_filt_study.csv")

# 3. try different subsets for model
