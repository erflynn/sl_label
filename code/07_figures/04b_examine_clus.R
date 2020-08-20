
#

require('tidyverse')
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv",
                          col_types="cccccccdld")
sample_source <- read_csv("data/sample_source_type.csv")
table(sample_source$source_type)
comb_metadata_w_src <- comb_metadata %>% 
  inner_join(sample_source %>% 
               select(acc, source_type), by=c("sample_acc"="acc"))

tissue_data <- comb_metadata_w_src %>%
  filter(source_type == "tissue" & !is.na(p_male))

stopifnot(nrow(comb_metadata_w_src)==nrow(comb_metadata))

ggplot(comb_metadata_w_src %>% 
         filter(source_type %in% c("named_cl", "unnamed_cl","tissue","primary_cells")) %>%
         mutate(source_type=fct_recode(source_type, "cell line"="named_cl", "primary cells"="primary_cells")) %>%
         mutate(source_type=fct_collapse(source_type, "cell line"=c("cell line", "unnamed_cl"))) %>%
         rename("sample source"=source_type), 
       aes(x=p_male, col=`sample source`)) +
  geom_density()+
  theme_bw()+
  facet_grid(organism~data_type, scales="free")+
  xlab("P(male)")
  
comb_metadata_w_src %>% group_by(organism, data_type, source_type) %>%
  count() %>%
  pivot_wider(names_from=source_type, values_from=n) %>%
  View()



tissue_data <- comb_metadata_w_src %>%
  filter(source_type == "tissue" & !is.na(p_male) & present & num_reads > 100000) %>%
  select(-platform, -source_type, -present, -num_reads)

tissue_data_lab <- tissue_data %>% filter(metadata_sex %in% c("male", "female"))

tissue_data %>% filter(metadata_sex == "mixed")

# ------ STOP ------ #

load("data/07_model_dat/human_microarray_sex_train_mat.RData")
# --> X_train, X_test, Y_train, Y_test
load("data/07_model_dat/fit_human_microarray_sex.RData") # --> fit
my_coef <- coef(fit, s=my.lambda)[,1]
nonzero_coef <- my_coef[my_coef!=0]
nz <- data.frame(nonzero_coef) %>% 
  rename(coef=nonzero_coef) %>% 
  arrange(coef) 
nz$gene <- rownames(nz)

f_genes <- nz %>% filter(coef < 0 & gene!="(Intercept)") %>% pull(gene) %>% head(5)
m_genes <- nz %>% filter(coef > 0 & gene!="(Intercept)") %>% pull(gene) %>% tail(5)

X_test2 <- X_test[,c(f_genes, m_genes)]
my_dist <- dist(X_test2)
preds <- predict(fit, newx=X_test, s=my.lambda, type="response")
df_pred <- data.frame(do.call(cbind,list("pred"=preds[,1], 
                 "acc"=rownames(preds), "true_lab"=Y_test))) 
df_pred %>%
  mutate(pred=as.numeric(as.character(pred))) %>%
  ggplot(aes(x=pred))+geom_histogram()+facet_grid(true_lab~.)
require('fpc')
pam.res <- pamk(my_dist, k=2, diss=TRUE)
table(pam.res$pamobject$clustering)


pcs <- prcomp(my_dist)
datf <- data.frame(cbind(pcs$rotation[,c("PC2", "PC3")], 
                         "cluster"=pam.res$pamobject$clustering),
                   "true_label"=Y_test)
ggplot(datf, aes(x=PC2, y=PC3))+
  geom_point(aes(shape=factor(cluster), color=Y_test), alpha=0.5)


# tsne

require('Rtsne')
tsne3 <- Rtsne(X_test, dims = 2, perplexity=15,theta=0.1, 
               is_distance=FALSE, verbose=TRUE,check_duplicates=FALSE, 
               max_iter = 1000)

tsne4 <- data.frame(cbind(tsne3$Y, "cluster"=Y_test))
colnames(tsne4) <-c("x1", "x2", "cluster")
tsne4$cluster <- factor(tsne4$cluster)
ggplot(tsne4, aes(x=x1, y=x2))+
  geom_point(aes(color=cluster))+ 
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
