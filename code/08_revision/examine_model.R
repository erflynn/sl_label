library(glmnet)
library(tidyverse)

# load models


# plot fit results

load_fold <- function(organism, data_type, idx){
  load(sprintf("data/06_fold_dat/fold_%s_%s_sex_%s.RData", organism, data_type, idx))
  tv_df <- tibble(do.call(rbind, lapply(res, function(x) x$tv)))
  y_df <- tibble(do.call(rbind, lapply(res, function(x) x$ydf)))
  tv_df$fold <- idx
  y_df$fold <- idx
  return(list("tv"=tv_df, "ydf"=y_df))
}
load_fold_dat <- function(organism, data_type){
  res <- lapply(1:6, function(i) load_fold(organism, data_type, i))
  tv_df <- tibble(do.call(rbind, lapply(res, function(x) x$tv)))
  y_df <- tibble(do.call(rbind, lapply(res, function(x) x$ydf)))
  tv_df$organism <- organism
  tv_df$data_type <- data_type
  y_df$organism <- organism
  y_df$data_type <- data_type
  return(list("tv"=tv_df, "ydf"=y_df))
}
hm_folds <- load_fold_dat("human", "microarray")
hr_folds <- load_fold_dat("human", "rnaseq")
mm_folds <- load_fold_dat("mouse", "microarray")
mr_folds <- load_fold_dat("mouse", "rnaseq")

plot_fold_breakdown <- function(df, organism, data_type){
  ggplot(df %>% unite("grp2", c(grp, alpha), remove=FALSE), 
         aes(x=factor(alpha), y=(1-class.1), col=grp, group=grp2))+
    geom_boxplot()+
    geom_point(position=position_dodge(0.75), alpha=0.5)+
    theme_bw()+
    xlab("elastic net parameter alpha")+
    ylab("accuracy")+
    ggtitle(sprintf("Hyperparameter tuning - %s %s", organism, data_type))
}

plot_fold_breakdown(hm_folds$tv, "Human", "microarray")
ggsave("figures/revision/hm_hyperparam.png")


# we noticed that where alpha=0

# We then computed the median classification error and 
# median lambda for each alpha across all six validation folds, 
# and selected the value of alpha (and its corresponding lambda)
# with the lowest median error.
# alpha = 1 --> lasso
# alpha = 0 --> ridge

grab_best <- function(df){
  df %>%
    group_by(alpha, grp) %>%
      summarize(median_acc=median(1-class.1),
                median_lambda=median(lambda)) %>%
      pivot_wider(names_from=grp, values_from=median_acc) %>%
      filter(alpha!=0) %>%
      ungroup() %>%
      filter(valid==max(valid))
}

hm_folds$tv %>% grab_best() 
mm_folds$tv %>% grab_best()
hr_folds$tv %>% grab_best()

# tie break - chose this one b/c it has a narrower range but it is pretty arbitrary
mr_folds$tv %>% grab_best() %>% filter(alpha==max(alpha))


train_model <- function(organism, data_type, best.alpha, best.lambda){

  if (data_type=="microarray"){
    standardizeFlag = TRUE
  } else {
    standardizeFlag = FALSE
  }
  
  # load the training data
  load(sprintf("data/07_model_dat/%s_%s_sex_train_mat.RData", organism, data_type))
  fit <- glmnet(X_train, Y_train, family="binomial", 
                alpha=best.alpha,
                standardize=standardizeFlag)
  
  #preds_train <- predict(fit, newx=X_train,  type="response")
  preds_class_train <- sapply(predict(fit, newx=X_train,  s=best.lambda, type="class"), as.numeric)
  sum(preds_class_train==Y_train)/length(Y_train) 
  
  preds_test <- predict(fit, newx=X_test,  s=best.lambda, type="response")
  preds_class_test <- sapply(predict(fit, newx=X_test, s=best.lambda, type="class"), as.numeric)
  sum(preds_class_test==Y_test)/length(Y_test) 

  my.lambda <- best.lambda; my.alpha <- best.alpha
  mat_coef <- coef(fit, s=my.lambda) %>% as.matrix()
  save(fit, my.alpha, my.lambda, 
       file=sprintf("data/07_model_dat/fit_%s_%s_%s.RData", organism, data_type, "sex"))
  nonzero_coef <- mat_coef[mat_coef[,1]!=0,]
  coef_df <- data.frame(cbind("gene"=names(nonzero_coef), coef=nonzero_coef)) 
  
  coef_df2 <- coef_df %>% arrange(as.numeric(as.character(coef)))
  return(coef_df2)
}

cf <- train_model("mouse", "microarray", mm_folds$tv %>% grab_best() %>% pull(alpha),
            mm_folds$tv %>% grab_best() %>% pull(median_lambda))
cf <- cf %>% mutate(coef=as.numeric(as.character(coef)),
                    gene=as.character(gene)) %>% 
  arrange(desc(coef)) %>% filter(gene !="(Intercept)")
cf %>% write_csv("data/07_model_dat/coef_df_mouse_microarray.csv")

mr_cf <- train_model("mouse", "rnaseq", mr_folds$tv %>% grab_best() %>% filter(alpha==max(alpha)) %>% pull(alpha),
                  mr_folds$tv %>% grab_best() %>% filter(alpha==max(alpha))%>% pull(median_lambda))
mr_cf <- mr_cf %>% mutate(coef=as.numeric(as.character(coef)),
                    gene=as.character(gene)) %>% 
  arrange(desc(coef)) %>% filter(gene !="(Intercept)")
mr_cf %>% write_csv("data/07_model_dat/coef_df_mouse_rnaseq.csv")



hm_cf <- train_model("human", "microarray", hm_folds$tv %>% grab_best() %>% pull(alpha),
                     hm_folds$tv %>% grab_best() %>% pull(median_lambda))
hm_cf <- hm_cf %>% mutate(coef=as.numeric(as.character(coef)),
                          gene=as.character(gene)) %>% 
  arrange(desc(coef)) %>% filter(gene !="(Intercept)")
hm_cf %>% write_csv("data/07_model_dat/coef_df_human_microarray.csv")

hr_cf <- train_model("human", "rnaseq", hr_folds$tv %>% grab_best() %>% pull(alpha),
                     hr_folds$tv %>% grab_best() %>% pull(median_lambda))
hr_cf <- hr_cf %>% mutate(coef=as.numeric(as.character(coef)),
                          gene=as.character(gene)) %>% 
  arrange(desc(coef)) %>% filter(gene !="(Intercept)")
hr_cf %>% write_csv("data/07_model_dat/coef_df_human_rnaseq.csv")



# --- TODO  --- # 
# try models with different types of genes
#  0) how many Y chr reads show up in females?
#    - stratify by PAR
#  1) all nonPAR
#  2) XCI + nonPAR Y
#  3) XCI + PAR X + all Y
#  4) what we did: everything
#  5) grab top genes (where are they?)
#  6) XIST, RPS4Y1, KDM5D




convert_genes <- getBM(attributes = c("external_gene_name",
                                      "ensembl_gene_id",
                                      "start_position",
                                      "end_position",
                                      "band",
                                      "chromosome_name"), 
                       filters = "ensembl_gene_id", values = cf %>% pull(gene),
                       mart=mart) %>%
  arrange(start_position)


cf_annot <- cf %>% left_join(convert_genes, by=c("gene"="ensembl_gene_id")) %>% 
  left_join(mouse_genes %>% dplyr::select(`Gene Symbol`, `Sex bias`, `Escapes XCI in`), 
            by=c("external_gene_name"="Gene Symbol"))
ggplot(cf_annot, aes(y=external_gene_name, x=coef, fill=`Sex bias`)) +
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("coefficient in model")+
  ylab("")

cf_annot$external_gene_name <- factor(cf_annot$external_gene_name, 
                                      levels=cf_annot$external_gene_name)
cf_annot2 <- cf_annot %>% mutate(`XCI status`=case_when(
  chromosome_name=="Y" ~ "Y chromosome",
  !is.na(`Escapes XCI in`) ~ "XCI escape",
  TRUE ~ "XCI"
))
ggplot(cf_annot2, aes(y=external_gene_name, x=coef, fill=`XCI status`)) +
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("coefficient in model")+
  ylab("")




examine_coef <- function(coef_df2){
  # load a gene map
  if (organism=="human"){
    load("../cersi_code/data/rb_ensembl_to_hgnc.RData")
  } else {
    # grab this!
    #load() 
  }
  
  if (data_type=="rnaseq"){
    # how do we map these?
    coef_df3 <- coef_df2 %>% 
      rename(transcript=gene) %>%
      left_join(xy_genes %>% select(-chromosome_name), by=c("transcript")) %>%
      left_join(ref_dat %>% select(all_of(ref_cols)) %>% unique(), 
                by=c("ensembl_gene_id"))
    
  } else {
    coef_df3 <- coef_df2 %>%
      left_join(convert_genes, #%>% select(all_of(ref_cols)) %>% unique(), 
                by=c("gene"="ensembl_gene_id")) %>%
      mutate(coef=as.numeric(as.character(coef)))
  }
  
  coef_df3 %>% write_csv(sprintf("data/07_model_dat/%s_%s_%s_coef.csv", prefix, data_type, ds))
}

# --- for HUMAN --- #
# add info from Tukianien abt PAR/escape genes

head(coef_df3 )
symbol_j <- coef_df3 %>%
  mutate(hgnc_symbol=ifelse(hgnc_symbol=="PUDP", "HDHD1", hgnc_symbol)) %>%
  filter(gene != "(Intercept)") %>% 
  left_join(supp13, by=c("hgnc_symbol"="Gene name")) %>%
  filter(!is.na(gene))
#id_j <- coef_df3 %>% 
#  filter(gene != "(Intercept)") %>% 
#  inner_join(supp13, by=c("gene"="gene_id2")) %>%
#  filter(!is.na(gene))
# ... these are the same

symbol_j %>% View()
symbol_j$`Reported XCI status`
symbol_j <- symbol_j %>%
  arrange(desc(coef))
symbol_j$hgnc_symbol <- factor(symbol_j$hgnc_symbol, 
                               levels=symbol_j$hgnc_symbol)
symbol_j2 <- symbol_j %>% 
  mutate(`XCI status`=case_when(
    chromosome_name=="Y" ~ "Y Chromosome",
    chromosome_name=="X" ~ paste("X", `Reported XCI status`, sep=" ")
  ))
ggplot(symbol_j2, aes(y=hgnc_symbol, x=coef, fill=`XCI status`)) +
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("coefficient in model")+
  ylab("")+
  scale_fill_manual(values=c("lightblue", "pink", "darkgray", 
                             "purple", "lightgreen"))
ggsave("figures/revision/hm_fit_coef.png")
# PLOT


# TODO - fix color scheme 
ggplot(symbol_j2 %>% 
         filter(chromosome_name!="Y") %>% 
         mutate(`Sex-bias in GTEx`=ifelse(`Sex-bias in GTEx`=="NA",
                                          "Unknown", `Sex-bias in GTEx`)), 
       aes(y=hgnc_symbol, x=coef, fill=`Sex-bias in GTEx`))+  
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("coefficient in model")+
  ylab("")
ggsave("figures/revision/hm_sexbias_gtex.png")

symbol_j2 %>% filter(`XCI across tissues`!="NA" | `XCI in single cells`!="NA") %>% 
  select(hgnc_symbol, coef, `XCI status`,
         `XCI across tissues`, 
         `XCI in single cells`)


# repeat for all four fits!

# COUNTS RNASEQ