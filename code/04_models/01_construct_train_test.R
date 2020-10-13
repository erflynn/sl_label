#
require('tidyverse')


# reorganize the training and testing data

# //TODO - update:
# INPUT:
# - `list_samples.csv` (organism, data_type, sample, study)
# - `sample_labels.csv` (sample, sex, cell line, src, drug, date)
# -` sample_qc.csv` (sample, study, present, rc)

comb_metadata <- read_csv("data/01_sample_lists/combined_human_mouse_meta_v2.csv",
                          col_types="cccccccdld")

# remove samples that fail QC 
labeled_meta <- comb_metadata %>% 
  filter(data_type=="microarray" | (present & num_reads > 100000)) %>% 
  select(sample_acc, organism, data_type, study_acc)



metadata_sex <- read_csv("data/01_sample_lists/mapped_sl_all.csv")
sample_src <- read_csv("data/updated_sample_src.csv")
cell_line <- read_csv("data/02_metadata/cell_line_mapping.csv")


# remove cell line data
labeled_meta_non_cl <- labeled_meta %>%
  filter(sample_src != "cell") 


# number of samples w labels
labeled_meta_non_cl %>% 
  group_by(organism, data_type, metadata_sex)  %>%
  count()

# //TODO: deal w study sharing


# number of studies with labels
labeled_meta %>% 
  group_by(organism, data_type)  %>%
  separate_rows(study_acc) %>%
  distinct(study_acc) %>%
  count()


# split 1/3 training, 2/3 testing
labeled_meta %>% 
  group_split()


# for data with metadata sex labels - divide into:
#  A train studies
#     - training samples
#     - extra train study 
#  B test studies
#     - test samples
#     - extra test study 
#  C extra?



# This leads to:
# 1. training - specific studies, but only use a subset of samples
# 2. testing - specific studies, only a subset of samples
# 3. extended test - specific studies 
#   - single sex
#   - mixed sex
#   - HC data
# 4. train studies
#   - train only samples
#   - all samples
# 5. test studies - all samples
# 6. *not in train studies <--  what we'll use for misannot detection*
# 7. anything not included?


# #### STOP #####
# 
# # --- update the RNA-seq 
# 
# my_files =list.files(pattern="fold_human_rnaseq_sex*")
# res_list <- do.call(rbind, lapply(my_files, function(my_f) {
#   load(my_f); 
#   res_df <- do.call(rbind, lapply(res, function(x) 
#     do.call(rbind, lapply(x, function(y) y$tv))))
#   return(res_df)
#   }))
# res_list %>% filter(grp=="valid") %>% arrange(class.1)
# 
# 
# res_list %>% filter(grp=="valid") %>% 
#   group_by(alpha, ngenes) %>% 
#   summarize(across(c(deviance.1:mae.1, lambda), median)) %>%
#   arrange(class.1)
# 
# # this actually looks good! esp for RNA-seq data. YAYYYY
# # --- for the cell line data not so much sadly --- #
# 
# # ---------------------------- #
# 
# # play with training and testing data...
# require('glmnet')
# 
# load("data/07_model_dat/human_rnaseq_sex_train_mat.RData")
# # --> X_test, X_train
# # --> Y_test, Y_train
# require('limma')
# 
# fit <- lmFit(t(X_train), as.numeric(Y_train))
# fit <- eBayes(fit)
# tt <-topTable(fit, number=ncol(X_train))
# tt$transcript <- rownames(tt)
# tt %>% head(NUM_GENES)
# # standardize == FALSE?