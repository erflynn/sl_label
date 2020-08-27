# 04d_mismatch_detect.R
# E Flynn
# 8/26/2020
# Code for mismatch detection
# 
# TODO:
# - sensitivity to parameters: scale, noise detection, cutoff

require('tidyverse')
require('mclust')

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", 
                          col_types="cccccccdcd")
by_study <- read_csv("data/study_sex_lab.csv")
source_type <- read_csv("data/sample_source_type.csv")

dat_src <- comb_metadata %>% 
  inner_join(source_type %>% select(acc, source_type), by=c("sample_acc"="acc"))
stopifnot(nrow(comb_metadata)==nrow(dat_src))

mixed_sex_studies <- by_study %>% 
  filter(label_type=="metadata" & 
           study_sex=="mixed sex" & 
           num_present >= 10)  

mixed_sex_samples <- dat_src %>% 
  separate_rows(study_acc, sep=";") %>%
  semi_join(mixed_sex_studies) %>%
  filter(!is.na(p_male) & metadata_sex!="unknown")

# --- function for extracting info from the mclust fit in a tibble ----- #
extractMcClustFit <- function(fit){
  G <- fit$G 
  mu <- fit$parameters$mean
  var <- sqrt(fit$parameters$variance$sigmasq)
  probs <- fit$z
  n <- nrow(probs)
  mu1 <- mu[[1]] ; var1 <- var[[1]]; p1 <- probs[,1]
  
  if (G==1){
    mu2=NA; p2=rep(NA, n); var2=NA
  } else {
    mu2 <- mu[[2]] ; var2 <- var[[2]]; p2 <- probs[,2]
  }
  tibble(class=fit$classification, p1=p1, p2=p2,
         G=rep(G, n), mu1=rep(mu1, n),
         mu2=rep(mu2, n), var1=rep(var1, n),
         var2=rep(var2, n))
}


splitEst <- function(dat){
  dat_split <- dat %>% 
    select(study_acc, organism, source_type,data_type, sample_acc, metadata_sex, p_male) %>% 
    group_split(study_acc) 
  my_fits <- lapply(dat_split, function(x) {
    predCol= x %>% select(p_male);
    nnoise=which(predCol>0.3 & predCol < 0.7);
    
    if (length(nnoise) < 1 | length(nnoise) > 0.33*nrow(x)){
      fit=Mclust(predCol, G=1:2, model="V",
                 prior=priorControl(scale=0.15));
    } else {
      fit=Mclust(predCol, G=1:2, model="V",
                 initialization= list(noise=nnoise),
                 prior=priorControl(scale=0.15));
    }
    df <- cbind(x, extractMcClustFit(fit))
    return(df)
  })
  my_fits2 <- do.call(rbind, my_fits) %>% as_tibble()
  return(my_fits2)
}

clusBreakdown <- function(x){
  x %>%
    distinct(study_acc, G) %>% 
    group_by(G) %>% 
    count() %>% 
    ungroup() %>% 
    mutate(frac=n/sum(n)) %>%
    rename(num_clus=G)
}


fits2 <- splitEst(mixed_sex_samples)

fits2 %>% clusBreakdown()# 26.9%

fit3 <- fits2 %>% 
  mutate(diff_mu=mu2-mu1) %>%
  mutate(study_sep=(mu1 < 0.5 & mu2 > 0.5 & diff_mu > 0.3)) %>%
  mutate(match=case_when(
    class==0 ~ "unclassified",
    (metadata_sex=="female" & p2 > 0.9) ~ "mismatch",
    (metadata_sex=="male" & p1 > 0.9) ~ "mismatch",
    (metadata_sex=="male" & p2 > 0.9) ~ "match",
    (metadata_sex=="female" & p1 > 0.9) ~ "match",
    TRUE ~ "unclear"))

fit4 <- fit3 %>% filter(G!=1 & study_sep) 
  
countFracMatch <- function(dat){
  dat %>%
    mutate(cell_line=(source_type %in% c("unnamed_cl", "named_cl"))) %>%
    group_by(organism, cell_line, data_type) %>%
    mutate(tot=n()) %>%
    group_by(organism, cell_line, data_type, match, tot) %>% 
    count() %>%
    pivot_wider(names_from=match, values_from=n, values_fill=0) %>%
    ungroup() %>%
    mutate(across(match:unclear, ~./tot)) %>%
    arrange(cell_line, organism, data_type)
}
fit4 %>% countFracMatch()

mismatch <- fit4 %>% filter(match=="mismatch") %>% distinct(study_acc)


# make a plot!!
fit4 %>% 
  filter(study_acc %in% (mismatch %>% sample_n(10) %>% pull(study_acc))) %>%
  ggplot(aes(x=p_male, y=metadata_sex, col=match))+
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
  facet_grid(study_acc~ .)+
  ylab("metadata sex")+
  xlab("P(male)")

ggsave("figures/paper_figs/new_mismatch_detect2.png")


# ----- look at single sex ----- #

f_studies <- by_study %>% 
  filter(label_type=="metadata" & study_sex=="female only" & num_present >= 8)   
m_studies <- by_study %>% 
  filter(label_type=="metadata" & study_sex=="male only" & num_present >= 8)  

m_samples <- dat_src %>% 
  separate_rows(study_acc, sep=";") %>%
  semi_join(m_studies) %>%
  filter(!is.na(p_male) & metadata_sex!="unknown")

f_samples <- dat_src %>% 
  separate_rows(study_acc, sep=";") %>%
  semi_join(f_studies) %>%
  filter(!is.na(p_male) & metadata_sex!="unknown")


matchBreakdown <- function(dat){
  dat %>%
    mutate(match=case_when(
      p_male < 0.3 & metadata_sex=="female" ~ "match",
      p_male > 0.7 & metadata_sex=="female" ~ "mismatch",
      p_male > 0.7 & metadata_sex=="male" ~ "match",
      p_male < 0.3 & metadata_sex=="male" ~ "mismatch",    
      TRUE ~ "unclear"
    )) %>%
    countFracMatch()
}

m_samples %>% matchBreakdown()
f_samples %>% matchBreakdown()
mixed_sex_samples %>% matchBreakdown()


# fit the models on the single sex data!
f_fits2 <- splitEst(f_samples)
m_fits2 <- splitEst(m_samples)
f_fits2 %>% clusBreakdown()
m_fits2 %>% clusBreakdown()

