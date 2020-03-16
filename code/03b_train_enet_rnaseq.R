
require('tidyverse')
require('glmnet')
require('data.table')

sl <- read_csv(sprintf("data/%s/02_sample_lists/%s_rnaseq_sex_lab.csv", 
                       prefix, prefix))

# then filter for what *is* present

set.seed(500)
if (nrow(sl) > 2000){
  sl_sm <- sl[sample(1:nrow(sl), 2000),]
} else {
  sl_sm <- sl  
}

# put together the data frame
grab_samples <- function(study_id){
  sample_list <- sl_sm %>% filter(study_acc==study_id)
  my.f <- sprintf("data/rnaseq/%s/01_study_mat/%s.csv", 
                  prefix, study_id)
  if (!file.exists(my.f) | file.info(my.f)$size==0){
    return(NA)
  }
  # what to do if the file doesn't exist
  dat <- fread(my.f, data.table=FALSE)
  if (ncol(dat) < 2){
    return(NA)
  }
  overlap.s <- intersect(colnames(dat), sample_list$sample_acc)
  #print(colnames(dat))
  if(!str_detect(colnames(dat)[1], "gene_name")){
    return(dat[,overlap.s])
  }
  colnames(dat)[1] <- c("gene_name")
  return(dat[,c("gene_name", overlap.s)])
 
}

#dat1 <- fread(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", 
#                     prefix, sl_sm$study_acc[[3]]), data.table=FALSE)
# if they're all the same length
sample_dat <- lapply(unique(sl_sm$study_acc), grab_samples)
sample_dat2 <- sample_dat[!is.na(sample_dat) & !is.null(sample_dat)]
row_lengths <- unlist(sapply(sample_dat2, nrow))
if (all(row_lengths==row_lengths[[1]])){
  study_d1 <- data.frame(do.call(cbind, lapply(sample_dat2, function(df) 
    data.frame(df) %>% select(-gene_name))))
  sample_expr_df <- cbind(sample_dat2[[1]]$gene_name, study_d1)
} else {
  print("Row lengths do not match" )
  sample_expr_df <- sample_dat2 %>% reduce(full_join, by="gene_name")
}

# if they're not all the same length
#rownames(sample_expr_df) <- dat1$gene_name
rownames(sample_expr_df) <- sample_expr_df$gene_name


# load the transcript to gene mapping
map_l <- fread(sprintf("data/rnaseq/transcriptome_indices/%s/long/genes_to_transcripts.txt", prefix), 
               header=FALSE, data.table = FALSE)
#map_s <- fread(sprintf("data/rnaseq/transcriptome_indices/%s/short/genes_to_transcripts.txt", prefix),
#               header=FALSE, data.table=FALSE)
## map_s is the same as map_l
colnames(map_l) <- c("gene", "transcript")


# add rownames
# make sure the data are ok
load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s_gene_map.RData", prefix))
xy_genes <- ref_dat %>% 
  filter(chromosome_name %in% c("X", "Y")) %>% 
  select(chromosome_name, ensembl_gene_id) %>% 
  unique()

xy_transcripts <- inner_join(xy_genes, map_l, by=c("ensembl_gene_id"="gene"))

gene_names <- rownames(sample_expr_df)
overlap.transcripts <- xy_transcripts %>% filter(transcript %in% gene_names)
write_csv(overlap.transcripts, sprintf("data/%s/04_sl_input/xy_genes_rnaseq.csv", prefix))

sample_expr_df2 <- sample_expr_df[overlap.transcripts$transcript,] %>% select(-gene_name)


# separate + shuffle
set.seed(310)
f_df_u <- sample_expr_df2[,sl_sm$sex=="female"]
f_df <- f_df_u[,sample(1:ncol(f_df_u), ncol(f_df_u)) ]
m_df_u <- sample_expr_df2[,sl_sm$sex=="male"]
m_df <-m_df_u[,sample(1:ncol(m_df_u), ncol(m_df_u)) ]
  
# train/test
expr_f_train <- f_df[,1:ceiling(ncol(f_df)/2)]
expr_f_test <- f_df[,(ceiling(ncol(f_df)/2)+1):ncol(f_df)]

expr_m_train <- m_df[,1:ceiling(ncol(m_df)/2)]
expr_m_test <- m_df[,(ceiling(ncol(m_df)/2)+1):ncol(m_df)]

expr_train <- cbind(expr_f_train, expr_m_train)
train_sex_lab <- c(rep(0, ncol(expr_f_train)), rep(1, ncol(expr_m_train)))

expr_test <- cbind(expr_f_test, expr_m_test)
test_sex_lab <- c(rep(0, ncol(expr_f_test)), rep(1, ncol(expr_m_test)))


# train model
# train lasso on these data
x_train <- as.matrix(t(expr_train))
x_train <- apply(x_train, c(1,2),as.numeric)

x_test <- as.matrix(t(expr_test))
x_test <- apply(x_test, c(1,2),as.numeric)

cvfit = cv.glmnet(x_train, train_sex_lab, family="binomial", alpha=0.5)
preds_train <- predict(cvfit, newx=x_train, s="lambda.1se", type="response")
preds_class_train <- sapply(predict(cvfit, newx=x_train, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_train==train_sex_lab)/length(train_sex_lab) 

preds_test <- predict(cvfit, newx=x_test, s="lambda.1se", type="response")
preds_class_test <- sapply(predict(cvfit, newx=x_test, s="lambda.1se", type="class"), as.numeric)
sum(preds_class_test==test_sex_lab)/length(test_sex_lab) 

mat_coef <- coef(cvfit, lambda="lambda.1se") %>% as.matrix()
nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef))
coef_df2 <- coef_df %>% left_join(xy_transcripts, by=c("gene"="transcript"))

coef_df2 %>% 
  arrange(coef) %>% 
  filter(!is.na(chromosome_name)) %>% 
  write_csv(sprintf("data/%s/04_sl_input/rnaseq_%s_coef.csv", prefix, prefix))

save(cvfit, file=sprintf("data/%s/04_sl_input/rnaseq_cvfit.RData", prefix))

