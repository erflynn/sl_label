
# cluster the distributions
head(study_grp)
res <- study_grp %>% 
  filter(!is.na(avg_p_study) & num_samples >= 3) %>%
  group_by(cl_acc) %>%
  mutate(n=n()) %>%
  filter(n >= 5) %>%
  group_split(cl_acc) 
res2 <- lapply(res, function(x) 
  data2hist(x$avg_p_study, algo="ManualBreaks", breaks=seq(0, 1, 0.05)))
names(res2) <- lapply(res, function(x) x$cl_acc[[1]])
mymat_study <- new("MatH", nrows=length(res2), ncols=1, 
             ListOfDist=res2, names.rows=names(res2), names.cols="density")
study_clusters <- lapply(2:10, function(i) WH_kmeans(mymat_study, k=i))
plot(x=2:10, y=lapply(study_clusters, function(x) x$solution$Crit)) # 6?
plot(x=2:10, y=lapply(study_clusters, function(x) x$quality)) # 6??

clus_df <- cbind("cl_acc"=sapply(res, function(x) x$cl_acc[[1]]), 
                 "cluster"=study_clusters[[5]]$solution$IDX) %>%
  as_tibble()
study_plt <- study_grp %>% 
  filter(cl_acc %in% clus_df$cl_acc) %>%
  left_join(clus_df) %>%
  left_join(df_switch %>% select(cl_acc, cl_annot_sex, cl_allele_sex, switching_category)) %>% 
  left_join(dip_test_col) %>%
  mutate(dip_pass=ifelse((dip_p > 0.05), "unimodal", "multimodal"))
head(study_plt)

study_plt %>% 
  group_by(cluster)  %>%
  summarize(num_f=sum(cl_allele_sex=="female"),
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
  ggplot(aes(x=avg_p_study, group=cl_acc, 
             col=`annotated switch`))+
  geom_density(alpha=0.8)+
  facet_grid(cluster~cl_annot_sex, scales="free")+
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


# 1. samples within studies

# cluster studies?
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


# 2. cell lines across all samples - #
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


# ---------- examples ----------- #
  
# https://stackoverflow.com/questions/40686753/cluster-histograms-using-earth-movers-distance-r
library(HistDAWass)
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

