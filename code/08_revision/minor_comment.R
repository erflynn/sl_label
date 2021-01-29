
library(tidyverse)

# grab the color scheme
library(RColorBrewer)
my.l <- brewer.pal(n = 8, name = 'Set2')
blues <- brewer.pal(9,name="Blues")
oranges <- brewer.pal(9,name="Oranges")
my.cols3 <- c (my.l[3],my.l[2], my.l[8])
my.cols4 <- c (my.l[3],my.l[4],my.l[2], my.l[8])
my.cols6 <- c (my.l[3], blues[4],my.l[4],oranges[4], my.l[2], my.l[8])

# load the data
full_ds <- read_csv("data/sample_metadata_filt.csv", col_types="cccccdldcc")




# ---- mostly m vs mostly f ---- #
study_breakdown <- full_ds %>% 
  separate_rows(study_acc, sep=";") %>%
  group_by(organism, data_type, label_type, study_acc) %>%
  summarize(n=n(),
         num_unlab=sum(sex_lab=="unlabeled"),
         num_male=sum(sex_lab=="male"),
         num_female=sum(sex_lab=="female")) %>%
  mutate(across(contains("num"), ~./n))
study_breakdown2 <- study_breakdown %>% 
  pivot_longer(c(num_male, num_female), names_to="sex", values_to="fraction") %>%
  mutate(sex=str_replace_all(sex, "num_", ""))

ggplot(study_breakdown2 %>% 
         filter(label_type=="metadata", num_unlab < 0.5 | 
                  (n >= 60 & n*(1-num_unlab) >= 30)),
       aes(x=fraction, col=sex))+
  geom_density()+
  facet_grid(data_type ~ organism)+
  theme_bw()+
  scale_color_manual(values=c(my.l[3], my.l[2]))+
  theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept=0.8, col="darkgray", lty=2)+
  geom_vline(xintercept=0.2, col="darkgray", lty=2)+
  xlab("study breakdown")+
  theme(legend.position = "None")
ggsave("figures/revision/study_sex_breakdown_metadata.png")


ggplot(study_breakdown2 %>% 
         filter(label_type=="expression", num_unlab < 0.5 | 
                  (n >= 60 & n*(1-num_unlab) >= 30)),
       aes(x=fraction, col=sex))+
  geom_density()+
  facet_grid(data_type ~ organism)+
  theme_bw()+
  scale_color_manual(values=c(my.l[3], my.l[2]))+
  theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept=0.8, col="darkgray", lty=2)+
  geom_vline(xintercept=0.2, col="darkgray", lty=2)+
  xlab("study breakdown")+
  theme(legend.position = "None")
ggsave("figures/revision/study_sex_breakdown_expression.png")


# fraction unlabeled plot
ggplot(study_breakdown %>%
        mutate(n=ifelse(n>500, 500, n)) %>%
         filter(label_type=="expression"),
       aes(x=n, y=num_unlab))+
  geom_point(alpha=0.3, position=position_jitter())+
  facet_grid(data_type ~ organism, scale="free")+
  theme_bw()+
  #scale_color_manual(values=c(my.l[3], my.l[2]))+
  theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept=60, col="darkgray", lty=2)+
  geom_hline(yintercept=0.5, col="darkgray", lty=2)+
  xlab("number of samples")+
  ylab("fraction unlabeled")
ggsave("figures/revision/unlabled_fraction_expression.png")

ggplot(study_breakdown %>%
         mutate(n=ifelse(n>500, 500, n)) %>%
         filter(label_type=="metadata"),
       aes(x=n, y=num_unlab))+
  geom_point(alpha=0.3, position=position_jitter())+
  facet_grid(data_type ~ organism, scale="free")+
  theme_bw()+
  #scale_color_manual(values=c(my.l[3], my.l[2]))+
  theme(panel.grid.minor = element_blank())+
  geom_vline(xintercept=60, col="darkgray", lty=2)+
  geom_hline(yintercept=0.5, col="darkgray", lty=2)+
  xlab("number of samples")+
  ylab("fraction unlabeled")
ggsave("figures/revision/unlabled_fraction_metadata.png")


# --- platform --- #


# ---- RNA-seq counts ---- #
head(full_ds)
rnaseq_ds <- full_ds %>% 
  filter(!is.na(num_reads), !is.na(present), data_type=="rnaseq" ) %>%
  select(-platform, -present, -data_type, -study_acc) %>%
  pivot_wider(names_from="label_type", values_from="sex_lab") 
ggplot(rnaseq_ds, aes(x=num_reads, fill=expression))+
  geom_histogram(alpha=0.9, bins=100)+
  facet_grid(organism~., scales="free")+
  theme_bw()+
  xlim(c(0, 40000000))+
  xlab("number of reads")+
  ylab("number of samples")
ggsave("figures/revision/reads_vs_label.png")

ggplot(rnaseq_ds, aes(x=num_reads, fill=(expression=="unlabeled")))+
  geom_histogram(position="fill", bins=100)+
  facet_grid(organism~., scales="free")+
  theme_bw()+
  theme(legend.position = "none")+
  xlim(c(0, 40000000))+
  ylab("fraction unlabeled (blue)")+
  xlab("number of reads")
ggsave("figures/revision/reads_vs_labeled_frac.png")
# FIX COLOR SCHEME
# TILE WINDOW
rnaseq_ds2 <- rnaseq_ds %>%
  mutate(is_labeled=(expression!="unlabeled"))

t.test(rnaseq_ds2$num_reads[rnaseq_ds2$is_labeled & rnaseq_ds2$num_reads < 300000], 
       rnaseq_ds2$num_reads[!rnaseq_ds2$is_labeled  & rnaseq_ds2$num_reads < 300000])

summary(lm(is_labeled ~ num_reads, data=rnaseq_ds2, family="binomial"))


# label similarity
rnaseq_ds %>% 
  filter(metadata!="unlabeled", expression != "unlabeled") %>%
  mutate(match=(expression==metadata)) %>%
  ggplot(aes(x=num_reads,fill=match))+ 
  geom_histogram(position="fill", bins=50)+
  theme_bw()+
  facet_grid(organism~., scales="free")+
  xlim(c(0, 60000000))+
  ylab("fraction of samples")+
  xlab("number of reads")
ggsave("figures/revision/reads_labels_match.png")

# fraction match

# fraction unlabeled
# -- this does NOT match
rnaseq_ds2 %>% 
  group_by(organism) %>% 
  mutate(n=n()) %>% 
  group_by(organism, expression, n) %>% 
  summarize(num=n(), frac=num/n) %>% 
  distinct()
# <60% of mouse labeled, <70% of human labeled

full_ds %>% 
  filter(!is.na(p_male) & data_type=="microarray") %>%
  filter(label_type=="expression") %>%
  group_by(organism) %>% 
  mutate(n=n()) %>% 
  group_by(organism, sex_lab, n) %>% 
  summarize(num=n(), frac=num/n) %>% 
  distinct()
# 87% of mouse, 90% of human labeled
# NEED TO ADD THIS TO THE OVERALL RESULTS from test datasets


# + number of each 

# some sort of correlation test    
  