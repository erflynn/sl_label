# grab all the sex labels

require('tidyverse')
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

dir.path <- sprintf("data/%s/05_sl_output/", prefix)
my.f <- list.files(path=dir.path)
res <- do.call(rbind, lapply(my.f, function(f) read_csv(sprintf("%s/%s", dir.path, f))))

res %>% write_csv(sprintf("data/%s/all_sl.csv", prefix))

metadata <- read.csv(sprintf("data/%s/02_sample_lists/%s_metadata.csv", prefix, prefix))
metadata2 <- metadata %>% 
  left_join(res, by=c("acc"="id")) %>%
  mutate(sex_lab=ifelse(pred < 0.5, "female", "male"))

metadata2 %>% filter(is.na(pred)) # this should be NONE
metadata2 %>% write.csv(file=sprintf("data/%s/02_sample_lists/%s_metadata2.csv", prefix, prefix), row.names=FALSE, quote=TRUE)

# 
# met_sex <- metadata2 %>% filter(!is.na(sex))
# met_sex_std <- met_sex %>% filter(sex %in% c("female", "male"))
# table(met_sex_std$sex, met_sex_std$sex_lab)
# 
# #female  male
# #female  50898  2684
# #male     4042 48334
# 
# sum(met_sex_std$sex==met_sex_std$sex_lab, na.rm=TRUE)/nrow(met_sex_std)
# # 0.94 in humans, 0.81 in mice
# met_sex_std2 <- met_sex_std %>% filter(!is.na(sex_lab)) # why are there a lot in mice
# sum(met_sex_std2$sex==met_sex_std2$sex_lab, na.rm=TRUE)/nrow(met_sex_std2)
# # 0.935 in mice
# 
#                                        
# met_sex2 %>% group_by(sex, sex_lab) %>% count() %>% arrange(desc(n))
# 
# met_sex2 <- met_sex %>% filter(!sex %in% c("female", "male"))
