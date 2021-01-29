# 01c_accuracy_platform.R
# E Flynn
# 08/11/2020
#
# Look at platform level accuracy
# Note - relies on extended test data generated in 01b_accuracy_counts.R
#
# Output:
#   - platform accuracy bar plot
#   - plot fraction labeled vs size / accuracy
#   - supplemental tables with breakdown of accuracy

require('tidyverse')
require('RColorBrewer')
blues <- brewer.pal(5,name="Blues")

sample_meta <- read_csv("data/01_sample_lists/rb_metadata/human_microarray_sample_metadata.csv",
                        col_types="cccccdcccccdd")


comb_metadata <- read_csv("data/sample_metadata_filt.csv", col_types="cccccdldcc")

# where does the platform data come from?
library('GEOmetadb')
con <- dbConnect(SQLite(), "../GEOmetadb.sqlite")
dbGetQuery(con, "SELECT gsm.gsm,gpl FROM gsm JOIN gse_gsm ON 
           gse_gsm.gsm=gsm.gsm WHERE gse IN ('GSE12039');")
dbGetQuery(con, "SELECT gse,gpl FROM gse_gpl WHERE gse IN  ('GSE12039');")

# TODO - have we removed the problem platforms?
#gse <- read_tsv("data/01_sample_lists/rb_metadata/gse12039.txt", col_names=F)
#gse[,c("X1", "X10")]

plat_list <- comb_metadata %>% 
  distinct(platform) %>%    
  mutate(full_str=platform) %>%
  mutate(platform=str_extract(platform, "(?<=\\().+(?=\\)$)")) %>%
  group_by(full_str) %>%
  mutate(platform=ifelse(
    str_detect(platform, "\\("), 
    str_split(platform, pattern="\\(")[[1]][[2]],
    platform)) %>%
  ungroup()
plat_list %>% group_by(platform) %>% 
  summarize(n=n(), full_str=paste(full_str, collapse=";")) %>% 
  filter(n>1) %>% View()

comb_metadata2 <- comb_metadata %>% 
  rename(full_str=platform) %>%
  left_join(plat_list, by=c("full_str"))

nrow(extended_test2) # 30559

platAcc <- function(my_organism, my_data_type){
  extended_test2 <- read_csv(sprintf("data/data_old/%s_%s_extended_test.csv", 
                                     my_organism, my_data_type))
  extended_test3 <- extended_test2 %>% left_join(comb_metadata2 %>% 
                                 distinct(study_acc, platform))
  df <- extended_test3 %>% filter(!is.na(metadata_sex), !is.na(expr_sex), 
                                   metadata_sex %in% c("male", "female"), 
                                   expr_sex %in% c("male", "female")) %>%
    mutate(unlab=(p_male > 0.3 & p_male < 0.7)) %>%
    mutate(match=(!unlab & metadata_sex==expr_sex)) %>%
    group_by(platform) %>%
    summarize(n=n(), accuracy=sum(match)/sum(!unlab), frac_lab=sum(!unlab)/n())
  df$organism <- my_organism
  df$data_type <- my_data_type
  return(df)
}



(hm_plat <- platAcc("human", "microarray"))
(hs_plat <- platAcc("human", "rnaseq"))
(mm_plat <- platAcc("mouse", "microarray"))
(ms_plat <- platAcc("mouse", "rnaseq"))

plat_dat <- do.call(rbind, list(hm_plat, hs_plat, mm_plat, ms_plat)) %>%
  arrange(organism, data_type, desc(accuracy)) 

plat_dat %>% filter(accuracy < 0.7)
plotPlatAcc <- function(ds, plot_title){
  ds %>%
    mutate(accuracy=ifelse(accuracy==0, 0.01, accuracy)) %>%
    arrange(desc(accuracy)) %>%
    #mutate(platform=str_extract(platform, "(?<=\\().+(?=\\))")) %>% # grab text in parens
    #group_by(platform) %>%
    #mutate(n2=1:n()) %>% # add a suffix if there are multiple w same name
    #ungroup() %>%
    #mutate(platform=ifelse(n2>1, paste(platform, n2, sep="-"), platform)) %>%
    mutate(platform=factor(platform, levels=unique(platform))) %>%
    mutate(num_samples=case_when(
      n < 10 ~ "<10",
      n < 50 ~ "10 - 49",
      n < 100 ~ "50 - 99",
      n < 500 ~ "100 - 499",
      TRUE ~ ">500",
    )) %>% 
    mutate(num_samples=factor(num_samples, levels=c("<10", "10 - 49", "50 - 99", "100 - 499", ">500"))) %>%
    ggplot(aes(x=platform, y=accuracy, fill=num_samples))+
    geom_bar(stat="identity")+
    theme_bw()+
    ggtitle(plot_title)+
    scale_fill_manual(values=blues)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    xlab("")
  
}

plotPlatAcc(hs_plat, "Human - RNA-seq")
ggsave("figures/paper_figs/supp_plat_acc_hr.png")

plotPlatAcc(hm_plat, "Human - Microarray")
ggsave("figures/paper_figs/supp_plat_acc_hm.png")

plotPlatAcc(ms_plat, "Mouse - RNA-seq")
ggsave("figures/paper_figs/supp_plat_acc_mr.png")

plotPlatAcc(mm_plat, "Mouse - Microarray")
ggsave("figures/paper_figs/supp_plat_acc_mm.png")

# write out the platform-level accuracy for a supplement
plat_dat <- do.call(rbind, 
                    list(hm_plat, hs_plat, mm_plat, ms_plat)) %>%
  arrange(organism, data_type, desc(accuracy))  %>%
  mutate(accuracy=ifelse(is.na(accuracy), 0, accuracy)) 

plat_dat %>% 
  select(organism, data_type, everything()) %>%
  rename(frac_labeled=frac_lab) %>%
  mutate(across(c(accuracy,frac_labeled), signif, 3)) %>%
  write_csv("tables/supp_plat_accuracy.csv")
plat_dat %>% distinct(organism, platform) #66
plat_dat %>% distinct(platform)
nrow(plat_dat)

plat_dat %>% 
  semi_join(plat_dat %>% filter(accuracy < 0.7), by="platform") 

######

# 13 were not included
missing_plat <- plat_list %>% 
  filter(!platform %in% plat_dat$platform) %>% 
  distinct(platform) %>%
  filter(!platform %in% 
           c("zebrafish", "drosophila2", "bovine", "rat2302", 
             "Illumina_RatRef-12_V1.0"))

comb_metadata2 %>% filter(platform %in% c("zebrafish", "drosophila2", "bovine", "rat2302", 
                                          "Illumina_RatRef-12_V1.0"))

# fraction of all samples tho!
plat_counts <- comb_metadata %>% 
  group_by(organism, data_type) %>% mutate(tot=n()) %>%
  group_by(organism, data_type, tot, platform) %>% count() %>%
  arrange(organism, data_type, desc(n)) %>%
  rename(num_samples=n) %>%
  mutate(frac_samples=num_samples/tot) %>%
  ungroup() %>%
  select(-tot)
plat_counts2 <- plat_counts %>% 
  left_join(plat_dat %>% 
              rename(num_in_test=n, frac_test_labeled=frac_lab, 
                     test_accuracy=accuracy), 
            by=c("organism", "data_type","platform")) 
plat_counts2 %>% filter(test_accuracy < 0.7)

plat_count_summary <- plat_counts2 %>%
  mutate(category=case_when(
    is.na(test_accuracy) ~ "not_in_test",
    (test_accuracy < 0.7) ~ "poor_accuracy",
    TRUE ~ "included"),
    approx_frac=ifelse(is.na(test_accuracy), 0, 
                       frac_samples*frac_test_labeled),
    approx_acc=frac_samples*test_accuracy) %>%
  select(-num_samples, -num_in_test) %>%
  group_by(organism, data_type, category) %>%
  summarize(frac_samples=sum(frac_samples),
            approx_frac=sum(approx_frac), 
            approx_acc=sum(approx_acc, na.rm=TRUE)) %>%
  pivot_longer(c(frac_samples, approx_frac, approx_acc), 
               names_to="stat", values_to="value") %>%
  filter((stat=="approx_frac" & category=="included")|
           (stat=="approx_acc" & category=="included") |
           !stat %in% c("approx_frac", "approx_acc")) %>%
  mutate(category=case_when(
    stat=="approx_frac" ~ "approx_frac",
    stat=="approx_acc" ~ "approx_acc",
    TRUE ~ category
  )) %>% 
  select(-stat) %>%
  unique() %>%
  pivot_wider(names_from=category, values_from=value, values_fill=0) 


plat_count_summary %>% 
  ungroup() %>%
  select(organism, data_type, poor_accuracy, not_in_test, 
         included, approx_acc, approx_frac) %>%
  rename(approx_accuracy=approx_acc,
         approx_frac_labeled=approx_frac) %>%
  mutate(across(where(is.numeric), signif, 4))  %>% 
  write_csv("tables/summary_platform.csv")

plat_dat %>%
  filter(frac_lab < 0.5)

# some have small n
plat_dat %>%
  filter(n <= 10)

ggplot(plat_dat, aes(x=n, y=accuracy))+
  xlim(0,500)+
  geom_point()+
  xlab("number of samples")+
  theme_bw()
ggplot(plat_counts2 %>%
         rename(`fraction of test set labeled`=frac_test_labeled) %>%
         mutate(across(everything(), ~replace_na(.,0))),
       aes(x=frac_samples, y=test_accuracy, 
           col=`fraction of test set labeled`))+
  #scale_x_continuous(trans=log_trans(),
  #                   breaks = trans_breaks("log", function(x) exp(x)),
  #                   labels = trans_format("log", function(x) signif(exp(x), 3)))+
  geom_point(alpha=0.7)+
  theme_bw()+
  xlab("fraction of all samples in platform")+
  ylab("accuracy in extended test set")
  
# LABEL THE ONES THAT WE'RE NOT USING

ggsave("figures/paper_figs/extra_plat_acc_by_size.png")
