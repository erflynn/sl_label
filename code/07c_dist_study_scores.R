load("data/summary_files.RData")
load("data/sample_level_sex.RData" )
h_summ %>% head()
h_map %>% head()
human2 %>% head()
h_summ2 <- h_summ %>% filter(num_samples >=8)


reformSamp <- function(dat){
  dat %>% select(study, num_samples, sex) %>%
    left_join(h_map, by=c("study"="study_acc")) %>% 
    left_join(human2, by=c("sample_acc"="acc")) %>%
    filter(!is.na(text_sex)) %>%
    rename(sample=sample_acc)
}
m_only <- h_summ2 %>% 
  filter(labeling_method=="metadata" & sex=="male-only") %>%
  reformSamp()
m_only_stat <- m_only %>% group_by(study) %>% 
  dplyr::summarize(mu=mean(pred), sigma=sd(pred))

ggplot(m_only %>% semi_join(m_only_stat %>% sample_n(10)), aes(x=pred))+
  geom_density()+
  facet_grid(rows=vars(study), scales="free")

f_only <- h_summ2 %>% filter(labeling_method=="metadata" & sex=="female-only") %>%
  reformSamp()
f_only_stat <- f_only %>% group_by(study) %>% 
  dplyr::summarize(mu=mean(pred), sigma=sd(pred))

ggplot(f_only %>% semi_join(f_only_stat %>% sample_n(10)), aes(x=pred, y=0))+
  geom_point(alpha=0.3, aes(col=factor(text_sex)))+
  facet_grid(rows=vars(study), scales="free")


mixed <- h_summ2 %>% 
  filter(labeling_method=="metadata" & sex=="mixed") %>%
  reformSamp()
mixed_stat <- mixed %>% group_by(study, text_sex) %>% 
  dplyr::summarize(mu=mean(pred), sigma=sd(pred)) %>%
  mutate(mu_l=mu-3*sigma, mu_u=mu+3*sigma) %>%
  mutate(mu_2l=mu-2*sigma, mu_2u=mu+2*sigma) %>%
  ungroup()
mixed_stat2 <- mixed_stat %>%
   pivot_longer(c("mu_2l", "mu_2u"), names_to="stat", values_to="stat_val") %>%
   mutate(stat=paste(stat,text_sex, sep="_"))  %>%
   select(-text_sex) %>%
   pivot_wider(id_cols = "study", names_from="stat", values_from="stat_val") %>%
   mutate(sex_sep=ifelse(mu_2u_female >= mu_2l_male, FALSE, TRUE))

likely_mislabeled <- 
  mixed %>% select(-num_samples, -sex) %>% 
  left_join(mixed_stat, by=c("study", "text_sex")) %>% 
  filter(pred < mu_2l | pred > mu_2u) %>%
  left_join(mixed_stat2 %>% select(study, sex_sep)) %>%
  filter(sex_sep) # 2152

# need to have separation between mu_f & mu_m
ggplot(mixed %>% semi_join(likely_mislabeled %>% select(study) %>% 
                             unique() %>% sample_n(10)),
                           aes(x=pred, y=0))+
         geom_point(alpha=0.3, aes(col=factor(text_sex)))+
         facet_grid(rows=vars(study), scales="free")


ggplot(mixed %>% semi_join(mixed_stat %>% sample_n(10)), aes(x=pred, y=0))+
  geom_point(alpha=0.3, aes(col=factor(text_sex)))+
  facet_grid(rows=vars(study), scales="free")



# write out the files we want!
m_only %>% select(sample) %>% unique() %>%  rename(acc=sample) %>%
  left_join(human_metadata %>% select(acc, idx, f_idx)) %>%
  arrange(idx) %>%
  write_csv("data/human_m_only_samp.csv")

f_only %>% select(sample) %>% unique() %>%  rename(acc=sample) %>%
  left_join(human_metadata %>% select(acc, idx, f_idx)) %>%
  arrange(idx) %>%
  write_csv("data/human_f_only_samp.csv")

mixed %>% select(sample) %>% unique() %>%  rename(acc=sample) %>%
  left_join(human_metadata %>% select(acc, idx, f_idx)) %>%
  arrange(idx) %>%
  write_csv("data/human_mixed_samp.csv")
