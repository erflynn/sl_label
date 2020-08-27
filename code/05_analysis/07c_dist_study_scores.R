load("data/summary_files.RData")
load("data/sample_level_sex.RData" )
h_summ %>% head()
human2 %>% head()
h_summ2 <- h_summ %>% filter(num_samples >=8)
h_map <- comb_metadata %>% filter(organism=="human") %>% select(study_acc, sample_acc) %>% separate_rows(study_acc, sep=";")

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
mixed_stat <- mixed %>% 
  group_by(study, text_sex) %>% 
  dplyr::summarize(mu=mean(pred), sigma=sd(pred)) %>%
  mutate(mu_l=mu-3*sigma, mu_u=mu+3*sigma) %>%
  mutate(mu_2l=mu-2*sigma, mu_2u=mu+2*sigma) %>%
  ungroup()


# ------- TRY MCLUST ------ #
require('mclust')


extractFit <- function(fit){
  G <- fit$G 

  mu <- fit$parameters$mean
  var <- fit$parameters$variance$sigmasq
  probs <- fit$z
  n <- nrow(probs)
  mu1 <- mu[[1]] ; var1 <- var[[1]]; p1 <- probs[,1]
  
  if (G==1){
    mu2=NA
    p2=rep(NA, n)
    var2=NA
  } else {
    mu2 <- mu[[2]] ; var2 <- var[[2]]; p2 <- probs[,2]
    
  }
  tibble(class=fit$classification, p1=p1, p2=p2,
         G=rep(G, n), mu1=rep(mu1, n),
         mu2=rep(mu2, n), var1=rep(var1, n),
         var2=rep(var2, n))
}

study_split <- mixed %>% 
  select(study, sample, text_sex, pred) %>% 
  group_split(study) 
fits <- lapply(study_split, function(x) 
  {predCol= x %>% select(pred);
  nnoise=which(predCol>0.4 & predCol < 0.6);
  if (length(nnoise) < 1){
    fit=Mclust(predCol, G=1:2, model="V",
               prior=priorControl(scale=0.15));
  } else {
    fit=Mclust(predCol, G=1:2, model="V",
               initialization= list(noise=nnoise),
               prior=priorControl(scale=0.15));
  }
  df <- cbind(x, extractFit(fit))
  return(df)})
fits2 <- do.call(rbind, fits) %>% as_tibble()

fit3 <- fits2 %>% 
  mutate(diff_mu=mu2-mu1) %>%
  mutate(study_sep=(mu1 < 0.5 & mu2 > 0.5 & diff_mu > 0.3)) %>%
  mutate(match=case_when(
         class==0 ~ "unclassified",
         (text_sex=="female" & p2 > 0.95) ~ "mismatch",
         (text_sex=="male" & p1 > 0.95) ~ "mismatch",
         (text_sex=="male" & p2 > 0.9) ~ "match",
         (text_sex=="female" & p1 > 0.9) ~ "match",
         TRUE ~ "unclear"))

nrow(fit3) # 77,966
fit4 <- fit3 %>% filter(G!=1 & study_sep) # 68,809
table(fit4$match)/nrow(fit4) # 2.70% mismatch, 2.60% unclass
mismatch <- fit4 %>% filter(match=="mismatch") %>% distinct(study)

"GSE30101"

fit4 %>% 
  #mutate(match_sex=paste(text_sex, match)) %>%
  filter(study %in% (mismatch %>% sample_n(10) %>% pull(study))) %>%
  ggplot(aes(x=pred, y=text_sex, col=match))+
  #geom_boxplot(aes(group=match_sex))+
  geom_point(alpha=0.3, position=position_jitter(0.01))+
  geom_vline(aes(xintercept=mu1), alpha=0.5)+
  geom_vline(aes(xintercept=mu2), alpha=0.8)+
  geom_vline(aes(xintercept=(mu2-1.96*var2)), col="gray", lty=2)+
  geom_vline(aes(xintercept=(mu2+1.96*var2)), col="gray", lty=2)+
  geom_vline(aes(xintercept=(mu1-1.96*var1)), col="gray", lty=2)+
  geom_vline(aes(xintercept=(mu1+1.96*var1)), col="gray", lty=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        strip.text.y.right = element_text(angle = 0))+
  scale_color_manual(values=c("dark blue", "red", "gray", "purple"))+
  facet_grid(study~ .)+
  ylab("metadata sex")+
  xlab("P(male)")
 
ggsave("figures/paper_figs/mismatch_detect6.png")
# stop

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
  filter(sex_sep) %>%
  filter(text_sex!=expr_sex) # 773

mis2 <-  mixed %>% select(-num_samples, -sex) %>% 
  left_join(mixed_stat, by=c("study", "text_sex")) %>%
  left_join(mixed_stat2 %>% select(study, sex_sep)) %>%
  filter(sex_sep) %>%
  filter(text_sex!=expr_sex)

mis3 <-  mixed %>% select(-num_samples, -sex) %>% 
  left_join(mixed_stat, by=c("study", "text_sex")) %>%
  left_join(mixed_stat2 %>% select(study, sex_sep)) %>%
  filter(!sex_sep) %>%
  filter(text_sex!=expr_sex)

likely_mislabeled2 <- 
  mixed %>% select(-num_samples, -sex) %>% 
  left_join(mixed_stat, by=c("study", "text_sex")) %>% 
  filter(pred < mu_l | pred > mu_u) %>%
  left_join(mixed_stat2 %>% select(study, sex_sep)) %>%
  filter(sex_sep) %>%
  filter(text_sex!=expr_sex) # 3673, # 4166

# need to have separation between mu_f & mu_m
ggplot(mixed %>% semi_join(likely_mislabeled %>% select(study) %>% 
                             unique() %>% sample_n(10)),
                           aes(x=pred, y=text_sex))+
         geom_point(alpha=0.3, aes(col=text_sex))+
         facet_grid(rows=vars(study), scales="free") +xlab("P(male)")+
  theme_bw()+ theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())


ggplot(mixed %>% semi_join(mis3 %>% select(study) %>% 
                             unique() %>% sample_n(10)),
       aes(x=pred, y=text_sex))+
  geom_point(alpha=0.3, aes(col=text_sex))+
  facet_grid(rows=vars(study), scales="free") +xlab("P(male)")+
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())


ggplot(tmp %>% semi_join(h_map %>% 
                           select(study_acc, sample_acc) %>% sample_n(10), 
                         by=c("acc"="sample_acc")),
       aes(x=pred, y=text_sex))+
  geom_point(alpha=0.3, aes(col=text_sex))+
  facet_grid(rows=vars(study_acc), scales="free") +xlab("P(male)")+
  theme_bw()+ theme( panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())


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


# filter to remove cell line
metadata <- read.csv("data/01_metadata/human_metadata.csv") 
no_cl <- metadata %>% filter(cl_line %in% c("", "--"))
h_map2 <- h_map %>% semi_join(no_cl , by=c("sample_acc"="acc"))
no_cl_h <- h_summ2 %>% filter(sex=="mixed") %>% 
  select(study, num_samples, sex) %>%
  inner_join(h_map2, by=c("study"="study_acc")) %>% 
  left_join(human2, by=c("sample_acc"="acc")) %>%
  filter(!is.na(text_sex)) %>%
  rename(sample=sample_acc)

text_pres <- no_cl_h %>% filter(text_sex %in% c("male", "female")) %>%
  unique()
text_pres %>% select(sample, expr_sex, text_sex) %>% 
  unique() %>% nrow() # 78557
text_pres %>% select(sample, expr_sex, text_sex) %>% 
  unique() %>% filter(expr_sex!=text_sex) %>% nrow() # 5858
# 7.5%


text_pres %>% select(sample, expr_sex, text_sex, pred) %>% 
  unique() %>% filter(expr_sex!=text_sex) %>% filter(pred > 0.7 | pred < 0.3) %>%
  nrow()
# 3709/78557 4.72%

text_pres %>% select(sample, expr_sex, text_sex, pred) %>% 
  unique() %>% filter(expr_sex!=text_sex) %>% filter(pred > 0.8 | pred < 0.2) %>%
  nrow()
# 2887/78557 = 3.67%

text_pres %>% select(sample, expr_sex, text_sex, pred) %>% 
  unique() %>% filter(expr_sex!=text_sex) %>% filter(pred > 0.9 | pred < 0.1)
# 1866 # = 2.38%



mis2_no <- mis2 %>% semi_join(no_cl, by=c("sample"="acc")) # 895
mis3_no <- mis3 %>% semi_join(no_cl, by=c("sample"="acc")) # 4227

my_dat <- text_pres %>% select(study, sample,text_sex, expr_sex, pred)

mixed_stat <- my_dat %>% group_by(study, text_sex) %>% 
  dplyr::summarize(mu=mean(pred), sigma=sd(pred)) %>%
  mutate(mu_l=mu-sigma, mu_u=mu+sigma) %>%
  mutate(mu_3l=mu-3*sigma, mu_3u=mu+3*sigma) %>%
  mutate(mu_2l=mu-2*sigma, mu_2u=mu+2*sigma) %>%
  ungroup()

my_dat2 <- my_dat %>% left_join(mixed_stat)

filt <- my_dat2 %>% 
  filter(expr_sex!=text_sex & pred > mu_3l & pred < mu_3u) # 4539

mixed_stat2 <- mixed_stat %>%
  pivot_longer(c("mu_l", "mu_u"), names_to="stat", values_to="stat_val") %>%
  mutate(stat=paste(stat,text_sex, sep="_"))  %>%
  select(-text_sex) %>%
  pivot_wider(id_cols = "study", names_from="stat", values_from="stat_val") %>%
  mutate(sex_sep=ifelse(mu_u_female >= mu_l_male, FALSE, TRUE))

filt %>% select(study, sample, pred) %>% left_join(mixed_stat2) %>% filter(sex_sep) %>%
  filter(pred < 0.1 | pred > 0.9)
# 2351


all_h <- human2 %>% semi_join(no_cl, by="acc") %>% unique()

n_tot <- all_h %>% filter(text_sex %in% c("male", "female")) %>% nrow()
n_m <- all_h %>% filter(text_sex %in% c("male", "female") & expr_sex != text_sex)  %>% nrow() 
nm1 <- all_h %>% filter(text_sex %in% c("male", "female") & expr_sex != text_sex )  %>% 
  filter(pred > 0.9 | pred < 0.1) %>% nrow() 
  
nm2 <- all_h %>% filter(text_sex %in% c("male", "female") & expr_sex != text_sex )  %>% 
    filter(pred > 0.8 | pred < 0.2) %>% nrow() 

tmp <- all_h %>% filter(text_sex %in% c("male", "female") & expr_sex != text_sex )  %>% 
  filter(pred > 0.8 | pred < 0.2)



## Testing for population probabilities
## Case A. Tabulated data
           # the same
x <- c(78557, 102261)
p <- c(1558, 1972)

chisq.test(x, p = p, rescale.p = TRUE)

