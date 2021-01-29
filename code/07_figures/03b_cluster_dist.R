# Code for clustering the distributions of cell line sex labels 

require('tidyverse')
require('HistDAWass')

# ----- PART ZERO ----- #
# define the dataset and cluster the study

annot_male <- df_switch %>% filter(cl_annot_sex=="male") %>% pull(cl_acc)

res0 <- study_grp %>% 
  filter(cl_acc %in% annot_male) %>%
  filter(!is.na(avg_p_study) & num_samples >= 3) %>%
  group_by(cl_acc) %>%
  mutate(n=n()) %>%
  filter(n >= 5) 

my_ds <- res0 %>% pull(cl_acc) # 87 cell lines
res <- res0%>%
  group_split(cl_acc) 

res2 <- lapply(res, function(x) 
  data2hist(x$avg_p_study, algo="ManualBreaks", breaks=seq(0, 1, 0.05)))
names(res2) <- lapply(res, function(x) x$cl_acc[[1]])
mymat_study <- new("MatH", nrows=length(res2), ncols=1, 
                   ListOfDist=res2, names.rows=names(res2), names.cols="density")
study_clusters <- lapply(2:10, function(i) WH_kmeans(mymat_study, k=i))


# QUALITY is much higher (0.92 vs 0.7) - is this significant?

plot(x=2:10, y=lapply(study_clusters, function(x) x$solution$Crit)) # 6?
plot(x=2:10, y=lapply(study_clusters, function(x) x$quality)) # 6??

clus_df <- cbind("cl_acc"=sapply(res, function(x) x$cl_acc[[1]]), 
                 "cluster"=study_clusters[[1]]$solution$IDX) %>%
  as_tibble()


study_plt <- study_grp %>% 
  filter(cl_acc %in% unique(my_ds)) #%>%
  #left_join(clus_df) %>%
  #left_join(df_switch %>% select(cl_acc, cl_annot_sex, cl_allele_sex, switching_category)) %>% 
  #left_join(dip_test_col) %>%
  #mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "multimodal"))
head(study_plt)

study_plt %>% 
  select(cl_acc, avg_p_study, se_p_study, cluster)


study_plt %>% 
  group_by(cluster)  %>%
  summarise(num_f=sum(cl_allele_sex=="female"),
            num_m=sum(cl_allele_sex=="male"),
            frac_f=num_f/(num_f+num_m),
            frac_uni=sum(dip_pass=="unimodal")/n())

length(unique(study_plt$cl_acc)) # 183
study_plt %>% 
  distinct(cl_acc, cluster, 
           cl_annot_sex, cl_allele_sex,
           switching_category) %>%
  filter(cl_annot_sex=="male") %>%
  group_by(cluster) %>% count()
# 5,6 = some switch
# 4 = hc f
# 2 = mostly f
# 3 = hc m
# 1 = mostly m
# 87 cell lines

study_plt %>% 
  mutate(`annotated switch`=(cl_allele_sex!=cl_annot_sex)) %>%
  mutate(cluster=factor(cluster, levels=c("4","2","5","6","1","3"))) %>%
  #filter(cl_annot_sex=="male") %>%
  #unite(cluster_sex, c(cluster, cl_annot_sex), remove=FALSE) %>%
  ggplot(aes(x=avg_p_study, group=cl_acc))+
  geom_density(alpha=0.3, col="gray")+
  facet_grid(cluster~cl_annot_sex, scales="free")+
  xlab("Avg Study P(male)")+
  #scale_color_manual(values=c("dark green", "purple"))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())

head(study_plt)
study_plt %>% group_by(cluster) %>%
  filter(!is.na(se_p_study)) %>%
  summarise(mean=mean(avg_p_study), sd=sd(avg_p_study))

table(clus_df$cluster)

clus_df %>% filter(cluster==1)

# cell_line, cell_name, allele_sex, number of refs, 
#  study_mean, study_sd
#  n_studies, n_samples
#  fraction of studies, fraction samples

# could also highlight multi-modal distribution?
# or draw violin plots?
study_plt2 <- study_plt %>%
  group_by(cl_acc) %>%
  mutate(mean_p=mean(avg_p_study, na.rm=T), 
         se_p=sd(avg_p_study),
         median_p=median(avg_p_study, na.rm=T), 
         num_studies=n(),
         
         median_nsamples=median(num_samples),
         mean_nsamples=mean(num_samples),
         sd_nsamples=sd(num_samples),
         num_samples=sum(num_samples)) %>%
  arrange(median_p) %>%
  left_join(cell_sex2 %>% select(cl_acc, cl_name)) 
study_plt2$cl_name <- 
  factor(study_plt2$cl_name, unique(study_plt2$cl_name))
study_plt2$cluster <- 
  factor(study_plt2$cluster, unique(study_plt2$cluster))




# ----- PART ONE ----- #
# ----- is the clustering better than random? ----- #

# random shuffle x 1000
ds <- study_grp %>% 
  filter(cl_acc %in% annot_male) %>%
  filter(!is.na(avg_p_study) & num_samples >= 3) %>%
  group_by(cl_acc) %>%
  mutate(n=n()) %>%
  filter(n>5) %>% 
  ungroup() %>% 
  select(-num_studies, -n, -study_acc) 

set.seed(12)

genRandomClus <- function(ds){
  null_res <- ds
  cl_vec <- null_res$cl_acc
  cl_vec_shuff <- cl_vec[sample(1:length(cl_vec),length(cl_vec))]
  null_res$cl_acc <- cl_vec_shuff
  null_res_split <- null_res %>% group_split(cl_acc)
  
  null_res2 <- lapply(null_res_split, function(x) 
    data2hist(x$avg_p_study, algo="ManualBreaks", 
              breaks=seq(0, 1, 0.05)))
  names(null_res2) <- lapply(null_res_split, function(x) 
    x$cl_acc[[1]])
  null_mat_study <- new("MatH", nrows=length(null_res2), ncols=1, 
                        ListOfDist=null_res2, 
                        names.rows=names(null_res2), 
                        names.cols="density")
  null_study_clusters = tryCatch({
    lapply(2:10, function(i) 
      WH_kmeans(null_mat_study, k=i))
  }, error = function(e){
    print(e)
    return(NA)
  })
  if(is.na(null_study_clusters)){
    return(NA)
  }
  # get criteria
  # get quality
  qual_vec <- sapply(null_study_clusters, function(x) x$quality)
  crit_vec <-  sapply(null_study_clusters, function(x) x$solution$Crit)
  #null_res_avg <- null_res %>% 
  #  group_by(cl_acc) %>% 
  #  summarise(avg_p2=mean(avg_p_study), se_p2=sd(avg_p_study))
  res_ds <- tibble("k"=2:10, "quality"=qual_vec, "crit"=crit_vec)
  return(res_ds)
}
set.seed(12)
null_res0 <- lapply(1:100, function(i) genRandomClus(ds))
null_df <- do.call(rbind, null_res0[!is.na(null_res0)])
head(null_df)
study_df <- tibble(k=2:10, 
       "quality"=sapply(study_clusters, function(x) x$quality),
       "crit"=sapply(study_clusters, function(x) x$solution$Crit)
)

# Quality of random vs what we have
ggplot(null_df, aes(x=k, y=quality))+
  geom_point(alpha=0.3)+
  geom_point(data=study_df, col="red")+
  ylim(0,1)+
  theme_bw()

ggplot(null_df, aes(x=factor(k), y=quality, group=k))+
  geom_violin()+
  stat_summary(aes(y = quality,group=k), 
               fun=mean, colour="black", geom="point")+
  geom_point(data=study_df, col="red")+
  ylim(0,1)+
  theme_bw()+
  xlab("number of clusters")+
  ylab("clustering quality (sum sq dev explained)")

ggsave("figures/paper_figs/clustering_quality.png")


# Generate random distributions
genRandomDist <- function(ds){
  null_res <- ds
  cl_vec <- null_res$cl_acc
  cl_vec_shuff <- cl_vec[sample(1:length(cl_vec),length(cl_vec))]
  null_res$cl_acc <- cl_vec_shuff
  null_res_avg <- null_res %>% 
    group_by(cl_acc) %>% 
    summarise(avg_p2=mean(avg_p_study), 
              se_p2=sd(avg_p_study))
    
  null_clus <- knn(train=res_avg[,c("avg_p2", "se_p2")], 
                     test=null_res_avg[,c("avg_p2", "se_p2")], 
                     k=3, cl=res_avg$cluster)
  return(list("dat"=null_res_avg, "clus"=null_clus))
}

rand_dist <- lapply(1:100, 
                    function(i) genRandomDist(ds))

list_dist <- do.call(rbind, lapply(1:100, function(x) 
  {df <- rand_dist[[x]]$dat; 
  df$idx <- rep(x, nrow(df));
  return(df)}))

# compare density average to null dist
ggplot(list_dist, 
       aes(x=avg_p2, group=idx))+
  geom_density(col="gray")+
  geom_density(data=res_avg %>% mutate(idx=1000), col="black")+
  theme_bw()+
  xlab("cell line average P(male) across studies")
ggsave("figures/paper_figs/cl_density_null.png")


# look at difference in counts per cluster
counts_per_clus <- 
  lapply(rand_dist, function(x) table(x$clus))
counts_per_clus
chisq.test(rbind(table(res_avg$cluster), 
                 counts_per_clus[[1]]))
chisq.test(rbind(counts_per_clus[[3]], 
                 counts_per_clus[[1]]))
df_tru <- data.frame(table(res_avg$cluster)) %>% 
  mutate(idx=1000) %>%
  rename(cluster=Var1, value=Freq)

data.frame(do.call(rbind, counts_per_clus)) %>%
  mutate(idx=1:length(counts_per_clus)) %>%
  pivot_longer(-idx, names_to="cluster") %>%
  mutate(cluster=str_replace(cluster, "X", "")) %>%
  ggplot(aes(x=cluster, y=value))+
  geom_boxplot()+
  geom_point(data=df_tru, col="red")+
  theme_bw()+
  ylab("number of cell lines")
ggsave("figures/paper_figs/null_clus_breakdown.png")



# 1. is the distribution different from random shuffling?
# ks.test() --> D statistic
# D statistic between each of the random

# 2. is the clustering quality better than we'd get from random?
#  quality for k=2:10 for each random shuffling

# 3. is the breakdown of samples to clusters different than random?
# knn --> predict cluster membership
# chisq test for categories



# objects needed:
# - study_grp
# - df_switch
# - dip_test_col
# - samples_cl_sl



#study_plt2$cl_acc <- 
#  factor(study_plt2$cl_acc, unique(study_plt2$cl_acc))


# ----- PART TWO ----- #
# ----- set up complicated cell line figure ----- #
# // TODO: make sure this all looks good

# A Overall distribution (randomly shuffled?)
# B-1A Fraction of samples
# [B-1B Fraction of studies]
# B-2 Boxplot studies
# B-3 Cluster
# B-4 Number of studies/samples
# [B-5 Dip-pvalue]
# [B-6 Reference]

# Idea:
#   Autosomal vs sex-related variation within the cell line?
#   tSNR for this

# A) overall distribution
study_plt2 %>% 
  filter(cl_acc %in% my_ds) %>%
  distinct(cl_name, mean_p, num_studies, num_samples) %>%
  ggplot(aes(x=mean_p))+geom_histogram(bins=20, alpha=0.7, col="black")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks=c(0,2,4,6,8,10))+
  ylab("Number of cell lines")+
  xlab("Cell Line Average Study P(male)")

require('gridExtra')

# B-1 Fraction of male samples

study_plt2 %>% distinct(cl_acc, cl_name) %>%
  filter(cl_acc %in% my_ds) %>%
  left_join(samples_cl_sl %>% filter(!is.na(cl_acc) & 
                                       cl_annot_sex=="male" &
                                       !is.na(p_male)) %>%
              select(sample_acc, study_acc, cl_acc, cl_name, expr_sex, p_male) %>%
              group_by(cl_acc) %>%
              summarise(num_samples=n(),
                        male=sum(p_male > 0.7),
                        female=sum(p_male < 0.3),
                        unlabeled=num_samples-male-female) %>%
              pivot_longer(cols=c(male, female, unlabeled), names_to="sex", values_to="count") %>%
              mutate(frac=count/num_samples)) %>%
  ggplot(aes(y=cl_name, x=frac, fill=sex))+
  geom_bar(stat="identity") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("fraction of samples")+
  ylab("cell name")+
  scale_fill_manual(values=my.cols3)

# B-2 study distribution

# //TODO: NAs are removed
study_plt2 %>% 
  filter(cl_acc %in% my_ds & !is.na(avg_p_study)) %>%
  #filter(cluster==4) %>%
  ggplot(aes(y=cl_name))+
  geom_boxplot(aes(x=avg_p_study),outlier.alpha=0.1)+
  #geom_violin(aes(y=avg_p_study), scale="count")+
  #geom_point(data=study_plt2 %>% filter(cluster==1) %>% 
  #             distinct(cl_name, mean_p, num_studies), 
  geom_point(data=study_plt2 %>% 
               distinct(cl_name, median_p, num_studies, num_samples, cluster),
             aes(x=median_p), 
             alpha=0.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #theme(axis.text.x = element_text(angle = 90, size=7,
  #                                 vjust = 0.5, hjust=1))+
  xlab("Average Study P(male)")+
  ylab("Cell line")+
  geom_vline(xintercept=0.5, alpha=0.5, lty=2)#+
  #facet_grid(.~cluster, scales="free")

# B-5 dip-test col
ggplot(study_plt2 %>% 
         distinct(cl_name, cl_acc, median_p, num_studies, num_samples, cluster) %>%
         left_join(dip_test_col), aes(y=cl_name, x=-log10(dip_p+0.000001)))+
  geom_point()+
  theme_bw()+
  geom_vline(xintercept = 
               -log10(0.05/length(unique(study_plt2$cl_acc))+0.000001), col="gray", lty=2, alpha=0.9)

# B-3 cluster
ggplot(study_plt2 %>% 
         distinct(cl_name, median_p, num_studies, num_samples, cluster),
       aes(y=cl_name, x=1, fill=cluster))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("")+
  ylab("cell line")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# B-6 Reference
allele_freq <- read_csv("data/00_db_data/cellosaurus_allele_freq.csv")

allele_freq2 <- allele_freq %>% 
  filter(`Not_detected`==0 & Y==0) %>%
  select(-Y, -Not_detected) %>%
  mutate(across(c(`X,Y`,`X`), ~./num_srcs)) %>%
  arrange(`X`) %>%
  mutate(cl_acc=tolower(cl_acc))

allele_freq3 <- allele_freq %>%
  select(-cl_name) %>%
  mutate(cl_acc=tolower(cl_acc)) %>%
  #mutate(across(c(`X,Y`,`X`, `Y`, `Not_detected`), ~./num_srcs)) %>%
  pivot_longer(cols=c(`X,Y`, `X`, `Y`, `Not_detected`), 
               names_to="alleles", values_to="count") 

# //TODO: NAs removed  
study_plt2 %>% 
  distinct(cl_acc, cl_name, median_p, num_studies, num_samples) %>%
  left_join(allele_freq3, by=c("cl_acc")) %>%
  ggplot(aes(y=cl_name, x=count, fill=alleles))+
  #aes(size=num_srcs))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("number of sources")+
  ylab("Cell line")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# B-4 Number of studies/samples

# coeff <- 10
# study_plt2 %>% 
#   distinct(cl_acc, cl_name, median_p, num_studies, num_samples) %>%
#   ggplot(aes(y=cl_name, x=num_samples))+
#   geom_bar(stat="identity")+
#   geom_point(aes(x=num_studies*coeff), col="red", alpha=0.5)+
#   ylab("Cell line")+
#   scale_x_continuous(
#     name="Number of samples",
#     sec.axis=sec_axis(~./coeff, name="Number of studies")
#   )+
#   theme_bw()

# --- make a df to write out --- #
study_plt2 %>% 
  distinct(cl_acc, cl_name, median_p, median_nsamples,
           sd_nsamples, num_studies, num_samples) %>%
  select(-median_nsamples, -sd_nsamples) %>%
  rename(median_p_male=median_p) %>%
  left_join(cl_level_lab %>% select(cl_acc, cl_allele_sex, num_f, num_m, num_unknown)) %>%
  select(cl_acc, cl_name, cl_allele_sex, everything()) %>%
  arrange(median_p_male) %>%
  mutate(median_p_male=round(median_p_male, 3)) %>%
  write_csv("tables/supp_cl_1018.csv")
  


study_plt2 %>% 
  distinct(cl_acc, cl_name, median_p, median_nsamples,
           sd_nsamples, num_studies, num_samples) %>%
  
  ggplot(aes(y=cl_name, x=num_studies))+
  geom_bar(stat="identity")+
  geom_point(aes(x=median_nsamples), col="red", alpha=0.5)+
  #geom_errorbarh(aes(xmin=mean_nsamples-1.96*sd_nsamples,
  #                   xmax=mean_nsamples+1.96*sd_nsamples), col="red", alpha=0.5)+
  ylab("Cell line")+
  xlab("Number of studies\nSamples per study (red)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



# ---- COMPARE TO CNVs in CCLE ---- #

cnv <- read_tsv("~/Downloads/CCLE_copynumber_byGene_2013-12-03.txt")

ccle <- read_tsv("~/Downloads/cl_annot") %>%
  mutate(Name=tolower(Name))

head(ccle$CCLE_ID)

cl_sex3 <- cell_sex2 %>% filter(cl_acc %in% unique(my_ds))

intersect(tolower(ccle$Name), tolower(cl_sex3$cl_name))
# 45!!
keep_cols <- ccle %>% filter(Name %in% tolower(cl_sex3$cl_name)) %>%
  pull(CCLE_ID)
head(cnv[,1:10])

cnv_df <- cnv[,c(c("EGID", "SYMBOL", "CHR"), keep_cols)]
cnv_df %>% 
  filter(CHR=="Y") %>% 
  select(-CHR) %>%
  pivot_longer(c(-EGID, -SYMBOL), 
               names_to="cell_line",values_to="expr") %>%
  ggplot(aes(x=expr, y=cell_line))+
  geom_boxplot()

cnv_df %>% 
  filter(!CHR %in% c("X", "Y")) %>% 
  select(-CHR) %>%
  pivot_longer(c(-EGID, -SYMBOL), 
               names_to="cell_line",
               values_to="expr") %>%
  ggplot(aes(x=expr, y=cell_line))+
  geom_boxplot(outlier.alpha=0.5)+
  theme_bw()+
  xlab("CNV")

cnv_df2 <- cnv_df %>% 
  filter(CHR %in% c("X", "Y")) %>%
  pivot_longer(c(-EGID, -SYMBOL, -CHR), 
               names_to="CCLE_ID",values_to="expr") %>%
  left_join(ccle %>% select(Name, CCLE_ID) %>%
              mutate(cl_name2=tolower(Name))) %>%
  right_join(study_plt2 %>% 
              mutate(cl_name2=tolower(cl_name)))


cnv_df2 %>%
  filter(is.na(CHR) | CHR=="X") %>%
  #filter(CHR %in% c("X", "Y")) %>%
  ggplot(aes(x=expr, y=cl_name))+
  geom_boxplot(outlier.alpha=0.1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("CCLE CNV")+
  #facet_wrap(~CHR)+
  geom_vline(xintercept=0,lty=2, alpha=0.5, col="gray")

comp_cnv <- cnv_df2 %>% 
  filter(!is.na(CHR)) %>%
  group_by(CHR, cl_name, median_p) %>% 
  summarise(median_cnv=median(expr)) 
comp_cnv_x <- comp_cnv %>% filter(CHR=="X")
comp_cnv_y <- comp_cnv %>% filter(CHR=="Y")

cor(comp_cnv_y$median_cnv, comp_cnv_y$median_p) 
# 0.774, pval=4.58 x 10**-10
cor(comp_cnv_x$median_cnv, comp_cnv_x$median_p) 
# -0.094, pval=0.537

comp_cnv %>%
  ggplot(aes(x=median_p, y=median_cnv))+
  geom_point(alpha=0.7)+
  xlab("Median CNV")+
  ylab("Median P(Male)")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  facet_wrap(~CHR)
plot(comp_cnv_x$median_cnv, comp_cnv_x$median_p) 
plot(comp_cnv_y$median_cnv, comp_cnv_y$median_p) 


cnv_df %>% 
  filter(CHR %in% c("X")) %>% 
  select(-CHR) %>%
  pivot_longer(c(-EGID, -SYMBOL), 
               names_to="CCLE_ID",values_to="expr") %>%
  left_join(ccle %>% select(Name, CCLE_ID) %>% mutate(cl_name2=tolower(Name))) %>%
  right_join(study_plt2 %>% 
               mutate(cl_name2=tolower(cl_name)))%>%
  ggplot(aes(x=expr, y=cl_name))+
  geom_boxplot(outlier.alpha=0.1)+
  theme_bw()+
  xlim(-8, 5)+
  xlab("X chromosome CNV")


f_cl <- ccle %>% filter(Gender=="female") %>% sample_n(50) %>% pull(CCLE_ID)
cnv_df2 <- cnv[,c(c("EGID", "SYMBOL", "CHR"), 
                  intersect(f_cl, colnames(cnv)))] 
cnv_df2 %>%
  filter(CHR %in% c("X")) %>% 
  select(-CHR) %>%
  pivot_longer(c(-EGID, -SYMBOL), 
               names_to="cell_line",values_to="expr") %>%
  ggplot(aes(x=expr, y=cell_line))+
  geom_boxplot(outlier.alpha=0.5)+
  theme_bw()+
  xlim(-8, 5)+
  xlab("CNV")


# -------------------------- STOP -------------------- #

# 2, 5, 1, 3


null_res_avg <- null_res %>% 
  group_by(cl_acc) %>% 
  summarise(avg_p2=mean(avg_p_study), se_p2=sd(avg_p_study))%>%
  left_join(null_clus_df)

res_avg <- res0 %>% 
  group_by(cl_acc) %>% 
  summarise(avg_p2=mean(avg_p_study), se_p2=sd(avg_p_study)) %>%
  left_join(clus_df)

library(class)

rand_clus_mod <- knn(train=res_avg[,c("avg_p2", "se_p2")], 
                test=null_res_avg[,c("avg_p2", "se_p2")], 
                k=3, cl=res_avg$cluster)
# // todo need to evaluate different
clus_mod <- knn(train=res_avg[,c("avg_p2", "se_p2")], 
                test=res_avg[,c("avg_p2", "se_p2")], 
                k=3, cl=res_avg$cluster)
table(rand_clus_mod)
table(clus_mod==res_avg$cluster)

pcs <- prcomp(res_avg[,c("avg_p2", "se_p2")])
clus_pc <- tibble("PC1"=pcs$x[,1], "PC2"=pcs$x[,2], "cluster"=res_avg$cluster)
ggplot(clus_pc, aes(x=PC1, y=PC2, col=cluster))+geom_point()+
  theme_bw()
ggsave("figures/paper_figs/clus_pcs.png")


pcs_rand <- prcomp(null_res_avg[,c("avg_p2", "se_p2")])
rand_clus_pc <- tibble("PC1"=pcs_rand$x[,1], "PC2"=pcs_rand$x[,2], 
                       "cluster"=null_res_avg$cluster)
ggplot(rand_clus_pc, aes(x=PC1, y=PC2, col=cluster))+geom_point()

pred_pc <- predict(pcs, null_res_avg[,c("avg_p2", "se_p2")])

ggplot(clus_pc, aes(x=PC1, y=PC2, col=cluster))+
  geom_point(alpha=0.7)+
  geom_point(data=data.frame(pred_pc), alpha=0.7, col="gray")+
  theme_bw()



ks.test(res_avg$avg_p2, null_res_avg$avg_p2) 
# different

plot(density(res_avg$avg_p2)) # multi-modal
plot(density(null_res_avg$avg_p2))
null_clus <- null_res %>% left_join(null_clus_df)
table(null_clus$cluster)
table(study_plt$cluster)



study_plt %>% filter(cluster==6)

ggplot( null_clus, 
       aes(x=avg_p_study, group=cl_acc))+
  geom_density(alpha=0.8)+
  facet_grid(cluster~., scales="free")+
  xlab("Avg Study P(male)")+
  scale_color_manual(values=c("dark green", "purple"))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())
       
# study_plt %>% filter(cluster==5) %>%
#   ggplot(aes(x=avg_p_study, group=cl_acc, col=dip_pass))+geom_density()+
#   facet_grid(cl_acc~., scales="free")+
#   theme_bw()+
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y=element_blank())

study_plt %>% filter(cluster==6) %>%
  ggplot(aes(x=avg_p_study, group=cl_acc, col=dip_pass))+geom_density()+
  facet_grid(cl_acc~., scales="free")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())


# -----  cluster distribution of samples within studies ------ #
my_s <- samples_cl_sl %>% 
  filter(!is.na(cl_acc) & !is.na(p_male)) %>%
  select(cl_acc, sample_acc, study_acc, p_male)  %>%
  separate_rows(study_acc, sep=";") %>%
  group_by(study_acc, cl_acc) %>%
  mutate(n=n()) %>% 
  filter(n>30) 

study_res <- my_s %>%
  ungroup() %>%
  group_split(study_acc) 
study_res2 <- lapply(study_res, function(x) 
  data2hist(x$p_male, algo="ManualBreaks", breaks=seq(0, 1, 0.05)))
names(study_res2) <- lapply(study_res, function(x) x$study_acc[[1]])
samples <- sample(1:length(study_res2),200)
study_res3 <- study_res2[samples]
names(study_res3) <- names(study_res2)[samples]

mymat4 <- new("MatH", nrows=length(study_res3), ncols=1, 
                   ListOfDist=study_res3, names.rows=names(study_res3), 
              names.cols="density")
study_clusters4 <- lapply(2:10, function(i) WH_kmeans(mymat4, k=i))

plot(x=2:10, y=lapply(study_clusters4, function(x) x$solution$Crit)) # 5?
plot(x=2:10, y=lapply(study_clusters4, function(x) x$quality)) # 5??

clus_df4 <- cbind("study_acc"=names(study_res3), 
                  "cluster"=study_clusters4[[5]]$solution$IDX) %>%
  as_tibble()

my_s4 <- my_s %>% 
  inner_join(clus_df4)

my_s4 %>%
  ggplot(aes(x=p_male, group=study_acc))+geom_density()+
  facet_grid(cluster~., scales="free")


my_s4 %>% filter(cluster==1) %>%
  ggplot(aes(x=p_male, group=study_acc))+geom_density()+
  facet_grid(study_acc~., scales="free")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())

my_s4 %>% filter(cluster==2) %>%
  ggplot(aes(x=p_male, group=study_acc))+geom_density()+
  facet_grid(study_acc~., scales="free")+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y=element_blank())

my_s4 %>% filter(cluster==4) %>%
  ggplot(aes(x=p_male, group=study_acc))+geom_density()+
  facet_grid(study_acc~., scales="free")


my_s4 %>% filter(cluster==5) %>%
  ggplot(aes(x=p_male, group=study_acc))+geom_density()+
  facet_grid(study_acc~., scales="free")




# --- cluster cell lines across all samples --- #
cl_samp <- samples_cl_sl %>% 
  filter(!is.na(cl_acc) & !is.na(p_male)) %>%
  select(cl_acc, sample_acc, p_male)  %>%
  group_by(cl_acc) %>%
  mutate(n=n()) %>% 
  filter(n>30) 
head(cl_samp)

samp_res <- cl_samp %>%
  ungroup() %>%
  group_split(cl_acc) 
samp_res2 <- lapply(samp_res, function(x) 
  data2hist(x$p_male, algo="ManualBreaks", breaks=seq(0, 1, 0.05)))
names(samp_res2) <- lapply(samp_res, function(x) x$cl_acc[[1]])
samples <- sample(1:length(samp_res2),200)
samp_res3 <- samp_res2[samples]
names(samp_res3) <- names(samp_res2)[samples]

mymat5 <- new("MatH", nrows=length(samp_res3), ncols=1, 
              ListOfDist=samp_res3, names.rows=names(samp_res3), 
              names.cols="density")
study_clusters5 <- lapply(2:10, function(i) WH_kmeans(mymat5, k=i))

plot(x=2:10, y=lapply(study_clusters5, function(x) x$solution$Crit)) # 7?
plot(x=2:10, y=lapply(study_clusters5, function(x) x$quality)) # 7??



clus_df5 <- cbind("cl_acc"=names(samp_res3), 
                  "cluster"=study_clusters5[[5]]$solution$IDX) %>%
  as_tibble()

my_s5 <- cl_samp %>% 
  inner_join(clus_df5)

my_s5 %>%
  ggplot(aes(x=p_male, group=cl_acc))+geom_density()+
  facet_grid(cluster~., scales="free")

# ---------------------------- STOP HERE ---------------- #

# ---------- examples ----------- #
  
# https://stackoverflow.com/questions/40686753/cluster-histograms-using-earth-movers-distance-r
v1 <- rnorm(n=100, mean = 10, sd = 1)  # cluster 1 (around 10)
v2 <- rnorm(n=100, mean = 50, sd = 5)  # cluster 2 (around 50)
v3 <- rnorm(n=100, mean = 100, sd = 10) # cluster 3 (around 100)
v4 <- rnorm(n=100, mean = 12, sd = 2)  # cluster 1
v5 <- rnorm(n=100, mean = 45, sd = 6)  # cluster 2
v6 <- rnorm(n=100, mean = 95, sd = 6)  # cluster 3

# label sample value 
d1_long <- do.call(cbind, list(v1, v2, v3, v4, v5, v6)) %>%
  as_tibble() %>%
  mutate(sample=1:n()) %>%
  pivot_longer(V1:V6, names_to="label", values_to="value") 

# create lists of histogram distributions
lod<-vector("list",6)
lod[[1]] <- data2hist(v1, type = "regular")
lod[[2]] <- data2hist(v2, type = "regular")
lod[[3]] <- data2hist(v3, type = "regular")
lod[[4]] <- data2hist(v4, type = "regular")
lod[[5]] <- data2hist(v5, type = "regular")
lod[[6]] <- data2hist(v6, type = "regular")

# combine separate lists into a matrix of histogram objects
mymat <- new("MatH", nrows=6, ncols=1, 
             ListOfDist=lod, names.rows=c(1:6), names.cols="density")

# calculate clusters pre-specifying number of clusters (k)
sol <- WH_kmeans(mymat, k=3)
clus_dat <- cbind(paste("V", names(sol$solution$IDX), sep=""), sol$solution$IDX) %>% 
  as_tibble() %>%
  rename(label=V1, cluster=V2)

d1_long %>% 
  left_join(clus_dat) %>%
  ggplot(aes(x=sample, y=value, col=cluster))+geom_point()+facet_grid(label~ .)


b1  <- c(rnorm(n=100, mean=9, sd=2) , rnorm(n=100, mean=200, sd=20))   # cluster 1 (around 10 and 200)
b2  <- c(rnorm(n=100, mean=50, sd=5), rnorm(n=100, mean=100, sd=10))  # cluster 2 (around 50 and 100)
b3  <- c(rnorm(n=100, mean=99, sd=8), rnorm(n=100, mean=175, sd=17)) # cluster 3 (around 100 and 175)
b4  <- c(rnorm(n=100, mean=12, sd=2), rnorm(n=100, mean=180, sd=40))  # cluster 1
b5  <- c(rnorm(n=100, mean=45, sd=6), rnorm(n=100, mean=80, sd=30))  # cluster 2
b6  <- c(rnorm(n=100, mean=95, sd=6), rnorm(n=100, mean=170, sd=25))  # cluster 3
b7  <- c(rnorm(n=100, mean=10, sd=1), rnorm(n=100, mean=210, sd=30))   # cluster 1 (around 10 and 200)
b8  <- c(rnorm(n=100, mean=55, sd=5), rnorm(n=100, mean=90, sd=15))  # cluster 2 (around 50 and 100)
b9  <- c(rnorm(n=100, mean=89, sd=9), rnorm(n=100, mean=165, sd=20)) # cluster 3 (around 100 and 175)
b10 <- c(rnorm(n=100, mean=8, sd=2), rnorm(n=100, mean=160, sd=30))  # cluster 1
b11 <- c(rnorm(n=100, mean=55, sd=6), rnorm(n=100, mean=110, sd=10))  # cluster 2
b12 <- c(rnorm(n=100, mean=105, sd=6), rnorm(n=100, mean=185, sd=21))  # cluster 3

d2_long <- do.call(cbind, list(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12)) %>%
  as_tibble() %>%
  mutate(sample=1:n()) %>%
  pivot_longer(V1:V6, names_to="label", values_to="value") 


lod<-vector("list",12)
lod[[1]] <- data2hist(b1, type = "regular")
lod[[2]] <- data2hist(b2, type = "regular")
lod[[3]] <- data2hist(b3, type = "regular")
lod[[4]] <- data2hist(b4, type = "regular")
lod[[5]] <- data2hist(b5, type = "regular")
lod[[6]] <- data2hist(b6, type = "regular")
lod[[7]] <- data2hist(b7, type = "regular")
lod[[8]] <- data2hist(b8, type = "regular")
lod[[9]] <- data2hist(b9, type = "regular")
lod[[10]] <- data2hist(b10, type = "regular")
lod[[11]] <- data2hist(b11, type = "regular")
lod[[12]] <- data2hist(b12, type = "regular")

mymat2 <- new("MatH", nrows=12, ncols=1, 
              ListOfDist=lod, names.rows=c(1:12), names.cols="density")

sol2 <- WH_kmeans(mymat2, k=6)

clus_dat2 <- cbind(paste("V", names(sol2$solution$IDX), sep=""), sol2$solution$IDX) %>% 
  as_tibble() %>%
  rename(label=V1, cluster=V2)

d2_long %>% 
  left_join(clus_dat2) %>%
  ggplot(aes(x=sample, y=value, col=cluster))+geom_point()+facet_grid(label~ .)

