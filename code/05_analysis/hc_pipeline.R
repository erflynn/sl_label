require('tidyverse')
require('bestNormalize')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- as.numeric(args[2])

MIN.READS <- 100000
# get the list of mixed sex studies
metadata_sex <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix),
                         stringsAsFactors = FALSE) %>% as_tibble() %>%
  select(acc, sex)

exp_sample <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))

sex_lab <- metadata_sex %>% filter(sex %in% c("male", "female"))
sex_lab2 <- exp_sample %>% left_join(sex_lab %>% rename(sample_acc=acc))
study_samp2 <- sex_lab2 %>% 
  mutate(qc_pass=(present & (num_reads >= MIN.READS))) %>%
  select(study_acc, sample_acc, sex, qc_pass) 

study_counts_sex <- study_samp2 %>%
  group_by(study_acc) %>% 
  summarize( num_f=sum(sex=="female", na.rm=TRUE), 
             num_m=sum(sex=="male", na.rm=TRUE), num_samples=n(), 
             num_qc=sum(qc_pass)) %>%
  mutate(study_sex=case_when(
    num_samples < 10 ~ "small", # less than 10 samples
    num_f==num_samples ~ "female-only",
    num_m==num_samples ~ "male-only", 
    (num_m > 5 & num_f > 5) ~ "mixed sex",
    (num_f + num_m == 0) ~ "unknown",
    TRUE ~ "other")) # not all of the samples pass

table(study_counts_sex$study_sex)
mixed_sex_studies <- study_counts_sex %>% filter(study_sex=="mixed sex")
mixed_sex_samples <- study_samp2 %>% semi_join(mixed_sex_studies) %>% unique()

# get the list of rows
load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s_gene_map.RData", prefix))
xy_genes <- read_csv(sprintf("data/rnaseq/%s/03_model_in/xy_genes_rnaseq.csv", prefix))

xy_genes2 <- xy_genes %>% 
  left_join(ref_dat %>% select(ensembl_gene_id, hgnc_symbol) %>% unique())
my_rows <- xy_genes2$transcript

# now load the studies in chunks
slStudy <- function(study_id){
  my_dat <- read_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))  %>% as.data.frame()
  rownames(my_dat) <- my_dat$gene_name
  my_dat$gene_name <- NULL 
  # select genes and rotate
  return(my_dat[intersect(my_rows,rownames(my_dat)),]) # return labels
}

# then boxcox transform
transformDat <- function(expr2){
  num_zeros <- apply(expr2, 1, function(x) sum(x==0)) # remove > 70% zeros
  expr3 <- expr2[(num_zeros/ncol(expr2) <= 0.7),]
  rowSums <- apply(expr3, 1, sum)
  boxcox_trans <- apply(expr3, 1, function(row) boxcox(row+0.5)$x.t)
  return(boxcox_trans)
}

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}
study.chunk <- extractChunk(mixed_sex_studies$study_acc, idx, 32)

expr2 <- do.call(cbind, lapply(study.chunk, slStudy))
expr4 <- transformDat(expr2)
expr4 <- t(expr4)

# get the genes for sex labeling
toker_genes <- xy_genes2 %>% 
  filter(hgnc_symbol %in% c("XIST", "KDM5D", "RPS4Y1"))  
massir_genes <- xy_genes2 %>% filter(chromosome_name=="Y")
toker_genes2 <- toker_genes %>% filter(transcript %in% rownames(expr4))
xist <- expr4[(toker_genes2 %>% filter(hgnc_symbol=="XIST"))$transcript,]
rps <- expr4[(toker_genes2 %>% filter(hgnc_symbol=="RPS4Y1"))$transcript,]
kdm <- expr4[(toker_genes2 %>% filter(hgnc_symbol=="KDM5D"))$transcript,]
massir_genes2 <- massir_genes %>% filter(transcript %in% rownames(expr4))

source("../sex_labeling/geo_pipeline/code/utils/sex_lab_utils.R")

# read in a study
studyAcc <- function(study_id){
  study_samples <- mixed_sex_samples %>% filter(study_acc==study_id)
  expr5 <- expr4[,study_samples$sample_acc]
  # cluster based on Toker
  
  macc <- massiRAcc(expr5, unlist(massir_genes))
  tacc <- tokerSexLab(expr5, f.genes=unlist(rownames(xist)), 
                      m.genes=c(unlist(rownames(rps)), unlist(rownames(kdm))))
  if (is.null(macc)) { macc <- rep(NA, nrow(study_samples))}
  
  if (is.null(tacc)) { tacc <- rep(NA, nrow(study_samples))}
  df <- data.frame(cbind(study_samples$sample_acc, study_samples$sex, macc, tacc))
  colnames(df) <- c("sample_acc", "metadata_sex", "massir_sex", "toker_sex")
  return(df)
}

study_res <- lapply(study.chunk, studyAcc)
study_df <- do.call(rbind, study_res)
study_df %>% write_csv(sprintf("data/rnaseq/%s/hc_%s.csv", prefix, idx))
exit()

study_df2 <- study_df %>% 
  mutate(matching=((massir_sex==metadata_sex & !is.na(massir_sex)) | 
                     (toker_sex==metadata_sex & !is.na(toker_sex))))
table(study_df2$matching)
study_df2 %>% filter(!matching)
# 12/215 do not match

# how do we compare?
imputed_sl <- read_csv(sprintf("data/02_labeled_data/%s_rnaseq_sl.csv", prefix)) %>%
mutate(imputed_sex=ifelse(pred > 0.5, "male", "female"))
study_df3 <- study_df2 %>% left_join(imputed_sl, by=c("sample_acc"="id"))
study_df3 %>% filter(!matching) %>% select(metadata_sex, massir_sex, imputed_sex)
hc <- study_df3 %>% filter(matching)
sum(hc$imputed_sex==hc$metadata_sex)/nrow(hc) # 100%...
