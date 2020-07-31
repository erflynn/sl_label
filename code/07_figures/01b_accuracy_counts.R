# 01b_accuracy_counts.R
# E Flynn
# 07/31/2020
#
# Code for calculating the accuracy (agreement) of our method
#
# TODO:
# - microarray - double check
# - get HC data and compare!
# - look at platform specificity
# - cutoffs

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

  acc=table(ds3$match)[["TRUE"]]/nrow(ds3)
  return(list("ds_assess"=ds_name, "cutoff"=cutoff, 
              "num_samples"=num_samples, "num_f"=num_f,
              "num_m"=num_m, "num_studies"=num_unique_studies, 
              "num_unlab_at_cutoff"=num_unlab,
              "accuracy"=acc))
}

createAccDf <- function(my_organism, my_data_type){
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
  print(length(intersect(train_studies, test_studies)) == 0)
  
  # filter to remove training samples (and studies)
  extended_test <- comb_metadata_filt %>%
    select(sample_acc, study_acc, metadata_sex, expr_sex, p_male) %>%
    anti_join(train) %>% # remove training samples
    group_by(sample_acc) %>%
    mutate(any_train=(length(intersect(str_split(study_acc, ";")[[1]], train_studies)) !=0)) %>%
    ungroup() 
  
  extended_test2 <- extended_test %>% filter(!any_train) 
  
  #
  # ---- set up for ms/ss ---- #
  by_study_sm <- by_study %>% 
    filter(label_type=="metadata" & !study_sex=="unknown" &
             organism==my_organism & data_type==my_data_type) %>%
    mutate(train_study=(study_acc %in% train_studies))
  
  
  
  # filter for a specific study size
  ms <- by_study_sm %>% filter(study_sex == "mixed sex" & num_tot >= 10)
  ss <- by_study_sm %>% filter(study_sex %in% c("male only", "female only") & num_tot >= 8)
  ms_sm <- by_study_sm %>% filter(study_sex == "mixed sex" & num_tot < 10)
  ss_sm <- by_study_sm %>% filter(study_sex %in% c("male only", "female only") & num_tot < 8)
  ms2 <- addSamples(ms, comb_metadata_filt)
  ss2 <- addSamples(ss, comb_metadata_filt)
  ms_sm2 <- addSamples(ms_sm, comb_metadata_filt)
  ss_sm2 <- addSamples(ss_sm, comb_metadata_filt)
  ms_v2 <- addSamples(ms %>% filter(!train_study), comb_metadata_filt)
  ss_v2 <- addSamples(ss %>% filter(!train_study), comb_metadata_filt)

  # // TODO add HC data 
  
  # calculate accuracies and put data frame together
  acc_df <- do.call(rbind, list(grabStats(train2, "train"),
                                grabStats(test2, "test"), 
                                grabStats(extended_test, "extended_test"),
                                grabStats(extended_test2, "extended_test_v2"),
                                grabStats(ss2, "single_sex"),
                                grabStats(ms2, "mixed_sex"),
                                grabStats(ss_sm2, "single_sex_sm"),
                                grabStats(ms_sm2, "mixed_sex_sm"),
                                grabStats(ss_v2, "single_sex_v2"),
                                grabStats(ms_v2, "mixed_sex_v2")))
  acc_df2 <- data.frame(acc_df)
  acc_df2$organism <- my_organism
  acc_df2$data_type <- my_data_type
  return(acc_df2)
}

hr_acc_df <- createAccDf("human", "rnaseq")
hm_acc_df <- createAccDf("human", "microarray")
mr_acc_df <- createAccDf("mouse", "rnaseq")
mm_acc_df <- createAccDf("mouse", "microarray")

acc_dat <- data.frame(do.call(rbind, list(hr_acc_df, hm_acc_df, mr_acc_df, mm_acc_df)))
acc_dat2 <- data.frame(apply(acc_dat , c(1,2), unlist))
acc_dat2 %>% select(organism, data_type, everything()) %>% write_csv("tables/s2_accuracies.csv")

# NOTE: microarray has data leakage!! eeps :/

# // TODO:
#  varying the threshold of p_male --> what is the accuracy
#  does this match the metric in that paper?
delt_cutoff <- c(0.5, 0.6, 0.7, 0.8, 0.9)

my_l <-do.call(rbind, list(grabStats(test2, "c0.5", 0.5),
     grabStats(test2, "c0.6", 0.6),
     grabStats(test2, "c0.7", 0.7),
     grabStats(test2, "c0.8", 0.8),
     grabStats(test2, "c0.9", 0.9)))
my_l

# ---- HC data ---- #
# only have it for microarray so far :/ 
#human_compare <- read_csv("data/human_sl_compare.csv")
#mouse_compare <- read_csv("data/mouse_sl_compare.csv")


# ---- platform specific? ---- #



# --- DEPRECATED --- #
# grab all the data we used for models! what was it and why?
# filter it out

# loadDat <- function(f,d) {
#   bind_rows(read_csv(sprintf("data/03_train/%s_%s_sex_pos_sample.csv", f, d)),
#             read_csv(sprintf("data/03_train/%s_%s_sex_neg_sample.csv", f, d))) 
# }
# 
# f.data <- lapply(c("human", "mouse"),
#        function(f){
#         bind_rows(
#           loadDat(f, "microarray") %>% rename(sample_acc=acc) %>% select(sample_acc),
#           loadDat(f, "rnaseq") %>% select(sample_acc)
#         )
#        }
# )
# train_data <- do.call(rbind, f.data)
# 
# # add in study, organism, data_type?
# train_data2 <- train_data %>% 
#   left_join(comb_metadata %>% select(sample_acc, study_acc, organism, data_type))
# 
# # be aware of multiple study data
# multi_study <- train_data2 %>% 
#   filter(str_detect(study_acc, ";")) %>%
#   distinct(study_acc) 
# 
# # ALL STUDIES - problem: I think we used too many studies!!!
# all_studies <- train_data2 %>% 
#   distinct(organism, data_type, study_acc) %>% 
#   separate_rows(study_acc, sep=";") %>%
#   unique()

