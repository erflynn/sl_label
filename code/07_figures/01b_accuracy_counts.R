# 01b_accuracy_counts.R
# E Flynn
# 07/31/2020
#
# Code for calculating the accuracy (agreement) of our method
# This also generates extended_test set files.
# 
# Produces:
# - evalulation datasets
# - accuracy supplement
# - accuracy at cutoff figure
# 
# TODO:
# - get HC data and compare!

require('tidyverse')
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")
by_study <- read_csv("data/study_sex_lab.csv")

# create a table w accuracies
#  organism | data_type | assess_ds | cutoff | num_samples | num_f | num_m | num_studies | accuracy

getUniqueStudies <- function(ds){
  ds %>% 
    distinct(study_acc) %>% 
    separate_rows(study_acc, sep=";") %>%
    unique() %>%
    pull(study_acc)
}

addSamples <- function(ds, comb_filt){
  comb_filt %>%
    separate_rows(study_acc, sep=";") %>%
    semi_join(ds) %>% 
    select(sample_acc, metadata_sex, expr_sex, p_male) %>% 
    unique() %>% 
    left_join(comb_filt %>% select(sample_acc, study_acc))
}


grabStats <- function(ds, ds_name, cutoff=0.5) {
  stopifnot(length(unique(ds$sample_acc))==nrow(ds))
  ds2 <- ds %>% filter(!is.na(metadata_sex), !is.na(expr_sex), 
                metadata_sex %in% c("male", "female"), 
                expr_sex %in% c("male", "female")) %>%
    mutate(match=(metadata_sex==expr_sex),
           cutoff_pass=
             ((expr_sex=="female" & (1-p_male) > cutoff) |
            (expr_sex=="male" & p_male > cutoff )))
  
  num_unlab <- ds2 %>% filter(!cutoff_pass) %>% nrow()
  ds3 <- ds2 %>% filter(cutoff_pass)
  num_samples <- nrow(ds3)
  num_unique_studies <- length(getUniqueStudies(ds3))
  ref_counts <- ds3 %>% group_by(metadata_sex) %>% count() 
  num_f <- ref_counts %>% filter(metadata_sex=="female") %>% pull(n)
  num_m <- ref_counts %>% filter(metadata_sex=="male") %>% pull(n)
  if (is.null(num_f)){ num_f=0 }
  if (is.null(num_m)){ num_m=0 }
  acc=(ds3 %>% filter(match) %>% count() %>%pull(n))/nrow(ds3)
  return(list("ds_assess"=ds_name, "cutoff"=cutoff, 
              "num_samples"=num_samples, "num_f"=num_f,
              "num_m"=num_m, "num_studies"=num_unique_studies, 
              "num_unlab_at_cutoff"=num_unlab,
              "accuracy"=acc))
}



createAccDf <- function(my_organism, my_data_type, cutoffs=NULL){
  load(sprintf("data/07_model_dat/%s_%s_sex_train_mat.RData", my_organism, my_data_type))
  # X_train, X_test, Y_train, Y_test
  
  train <- data.frame(cbind("sample_acc"=rownames(X_train), Y_train)) %>% as_tibble()
  test <- data.frame(cbind("sample_acc"=rownames(X_test), Y_test)) %>% as_tibble()
  
  comb_metadata_filt <- comb_metadata %>% 
    filter(organism==my_organism & data_type==my_data_type)
  
  train2 <- train %>% inner_join(comb_metadata_filt)
  test2 <- test %>% inner_join(comb_metadata_filt) 
  
  # grab the studies 
  train_studies <- train2 %>% getUniqueStudies()
  test_studies <- test2 %>% getUniqueStudies()
  print(length(intersect(train_studies, test_studies)))
  
  # filter to remove training samples (and studies)
  extended_test <- comb_metadata_filt %>%
    select(sample_acc, study_acc, metadata_sex, expr_sex, p_male) %>%
    anti_join(train) %>% # remove training samples
    group_by(sample_acc) %>%
    mutate(any_train=(length(intersect(str_split(study_acc, ";")[[1]], 
                                       train_studies)) !=0)) %>%
    ungroup() 
  
  extended_test2 <- extended_test %>% filter(!any_train) 
  write_csv(extended_test2, sprintf("data/%s_%s_extended_test.csv", my_organism, my_data_type))
  
  # ---- set up for ms/ss ---- #
  by_study_sm <- by_study %>% 
    filter(label_type=="metadata" & !study_sex=="unknown" &
             organism==my_organism & data_type==my_data_type) %>%
    mutate(train_study=(study_acc %in% train_studies))
  
  # filter for a specific study size
  ms <- by_study_sm %>% 
    filter(study_sex == "mixed sex" & num_tot >= 10) %>%
    filter(!train_study)
  
  ss_f <- by_study_sm %>% 
    filter(num_tot==num_f & num_tot >= 8)  %>%
    filter(!train_study)
  
  ss_m <- by_study_sm %>% 
    filter(num_tot==num_m & num_tot >= 8)  %>%
    filter(!train_study)
  ms2 <- addSamples(ms, comb_metadata_filt)
  ss_f2 <- addSamples(ss_f, comb_metadata_filt)
  ss_m2 <- addSamples(ss_m, comb_metadata_filt)

  # // TODO add HC data 

  if (is.null(cutoffs)){
    cutoffs <- c(0.5)
  }
  c_acc_dfs <- lapply(cutoffs, function(cutoff){
      c_acc_df <- do.call(rbind, list(grabStats(train2, "train", cutoff),
                                    grabStats(test2, "test", cutoff),
                                    grabStats(extended_test2, "extended_test", cutoff),
                                    grabStats(ms2, "mixed_sex", cutoff),
                                    grabStats(ss_f2, "single_sex_f", cutoff),
                                    grabStats(ss_m2, "single_sex_m", cutoff)
                                    ))
      return(c_acc_df)
    })
  acc_df <- do.call(rbind, c_acc_dfs)
  
  acc_df2 <- data.frame(acc_df)
  acc_df2$organism <- my_organism
  acc_df2$data_type <- my_data_type
  return(acc_df2)
}


hr_acc_df <- createAccDf("human", "rnaseq", c(0.5, 0.6, 0.7, 0.8, 0.9))
hm_acc_df <- createAccDf("human", "microarray", c(0.5, 0.6, 0.7, 0.8, 0.9))
mr_acc_df <- createAccDf("mouse", "rnaseq", c(0.5, 0.6, 0.7, 0.8, 0.9))
mm_acc_df <- createAccDf("mouse", "microarray", c(0.5, 0.6, 0.7, 0.8, 0.9))

lapply(list(hr_acc_df, hm_acc_df, mr_acc_df, mm_acc_df), function(x) apply(x, 1, unlist))
acc_dat <- data.frame(do.call(rbind, list(hr_acc_df, hm_acc_df, mr_acc_df, mm_acc_df)))
acc_dat[sapply(acc_dat[,"num_f"], function(x) length(x)==0),"num_f"]<-0
acc_dat[sapply(acc_dat[,"num_m"], function(x) length(x)==0),"num_m"]<-0
acc_dat2 <- data.frame(t(apply(acc_dat, 1, unlist)))


acc_dat3 <- acc_dat2 %>% select(organism, data_type, everything()) %>% arrange(cutoff, organism, data_type, ds_assess) 
acc_dat2 %>% 
  select(organism, data_type, everything()) %>%
  filter(cutoff==0.5) %>%
  select(organism, data_type, ds_assess, num_samples, num_f, num_m, num_studies) %>%
  distinct() %>%
  write_csv("tables/s2a_eval_ds.csv")

acc_dat3 %>%
  select(-num_studies, -num_f, -num_m) %>%
  as_tibble() %>%
  mutate(across(cutoff:accuracy, ~as.numeric(as.character(.)))) %>%
  mutate(num_total=num_samples+num_unlab_at_cutoff) %>%
  mutate(frac_labeled=num_samples/num_total) %>%
  select(-num_samples, -num_unlab_at_cutoff, -num_total) %>%
  mutate(across(c(accuracy,frac_labeled), signif, 3)) %>%
  write_csv("tables/s2b_accuracies.csv")


acc_dat3 <- acc_dat2 %>% as_tibble() %>%
  mutate(across(c(-ds_assess, -organism, -data_type), ~as.numeric(as.character(.)))) %>%
  mutate(frac_labeled=num_samples/(num_samples+num_unlab_at_cutoff),
         threshold_score=factor(cutoff))

cut <- acc_dat3 %>% 
  filter(cutoff==0.7) %>% 
  select(organism, data_type, ds_assess, accuracy, frac_labeled) %>%
  unite(organism_d,c(organism, data_type), sep=" - ") %>%
  mutate(across(c(accuracy,frac_labeled), signif, 3))

cut %>%  
  pivot_wider(id_cols=ds_assess, values_from=accuracy, names_from=organism_d) %>%
  write_csv("data/cutoff_accuracy.csv")

cut %>%
  pivot_wider(id_cols=ds_assess, values_from=frac_labeled, names_from=organism_d) %>% 
  write_csv("data/cutoff_labeled.csv")

acc_dat3 %>%
  filter(ds_assess != "train") %>%
  mutate(ds_assess=fct_recode(ds_assess,
    "single_sex (f)" = "single_sex_f",
    "single_sex (m)" = "single_sex_m"
  )) %>%
  mutate(ds_assess=fct_relevel(
    ds_assess, "test"
  )) %>% 
  rename(dataset=ds_assess) %>%
  ggplot(aes(x=frac_labeled, y=accuracy, col=threshold_score, group=threshold_score))+
  xlab("fraction of data labeled")+
  ylim(0.5, 1)+
  theme_bw() + 
  geom_vline(xintercept = 0.7, col="gray")+
  geom_vline(xintercept = 0.8, col="gray")+
  geom_vline(xintercept = 0.9, col="gray")+
  geom_hline(yintercept = 0.98, col="gray")+
  geom_hline(yintercept = 0.95, col="gray")+
  geom_hline(yintercept = 0.9, col="gray")+
  geom_point(alpha=0.8, aes(shape=dataset))+ 
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  facet_grid(data_type~organism)
# make a pretty plot
ggsave("figures/paper_figs/supp_accuracy_p_cutoff.png")


# ---- HC data ---- #
# only have it for microarray so far :/ 
#human_compare <- read_csv("data/human_sl_compare.csv")
#mouse_compare <- read_csv("data/mouse_sl_compare.csv")

