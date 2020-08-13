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


platAcc <- function(my_organism, my_data_type){
  extended_test2 <- read_csv(sprintf("data/%s_%s_extended_test.csv", 
                                     my_organism, my_data_type))
  extended_test3 <- extended_test2 %>% left_join(comb_metadata %>% 
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



plotPlatAcc <- function(ds, plot_title){
  ds %>%
    arrange(desc(accuracy)) %>%
    mutate(platform=str_extract(platform, "(?<=\\().+(?=\\))")) %>% # grab text in parens
    group_by(platform) %>%
    mutate(n2=1:n()) %>% # add a suffix if there are multiple w same name
    ungroup() %>%
    mutate(platform=ifelse(n2>1, paste(platform, n2, sep="-"), platform)) %>%
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
plat_dat <- do.call(rbind, list(hm_plat, hs_plat, mm_plat, ms_plat)) %>%
  arrange(organism, data_type, desc(accuracy))  %>%
  mutate(accuracy=ifelse(is.na(accuracy), 0, accuracy)) 

plat_dat %>% 
  select(organism, data_type, everything()) %>%
  rename(frac_labeled=frac_lab) %>%
  mutate(across(c(accuracy,frac_labeled), signif, 3)) %>%
  write_csv("tables/supp_plat_accuracy.csv")


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


not_test <- plat_counts2 %>% filter(is.na(test_accuracy)) %>%
  group_by(organism, data_type) %>%
  summarize(frac_not_tested=sum(frac_samples))

frac_samples <- plat_counts2 %>% filter(!is.na(test_accuracy) & test_accuracy > 0.7) %>%
  group_by(organism, data_type) %>%
  summarize(frac_samples_included=sum(frac_samples))

lab <- plat_counts2 %>% 
  filter(!is.na(frac_test_labeled) & test_accuracy > 0.7) %>%
  mutate(approx_frac=frac_samples*frac_test_labeled) %>%
  group_by(organism, data_type) %>%
  summarize(approx_frac_labeled=sum(approx_frac))

plat_res <- plat_counts2 %>% filter(test_accuracy < 0.7) %>%
  group_by(organism, data_type) %>%
  summarize(frac_poor_acc=sum(frac_samples)) %>%
  full_join(not_test) %>%
  full_join(frac_samples) %>%
  full_join(lab) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), signif, 4)) 
  
plat_res%>% write_csv("tables/summary_platform.csv")

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
  
ggsave("figures/paper_figs/extra_plat_acc_by_size.png")
