# 05b_drug_diff.R
# E Flynn
#
# Exploratory look at sex differences in drug response in studies that are sufficiently powered
#
# Tried NANOVA (package is horrible) and robust ANOVA. Got them both working.
# Plotted results nicely! So useful if we find anything lol. 
#
# TODO:
# - look at with positive control
# - estimate power
# - try ARTool <-- Gao et al paper suggests this is better
# - try discovery/validation split

#install.packages("~/Downloads/TANOVA_1.0.0.tar.gz", repos=NULL)
#require('TANOVA') # BAD
require('GEOquery')
require('limma')
require('tidyverse')
require('fdrtool')
require('Rfit') # for rank anova
require('biomaRt')


# ---- are there any studies where we could look at sex diff in drug response? ---- $
exp_metadata <- read_csv("data/01_metadata/human_experiment_metadata.csv", 
                         col_types="cccc") %>%
  bind_rows(read_csv("data/01_metadata/mouse_experiment_metadata.csv", 
                         col_types="cccc")) %>%
  bind_rows(read_csv("data/01_metadata/human_rnaseq_experiment_metadata.csv", 
                     col_types="cccc")) %>%
  bind_rows(read_csv("data/01_metadata/mouse_rnaseq_experiment_metadata.csv", 
                     col_types="cccc"))
study_db2 <- read_csv("data/drug_studies.csv")
large_drug_studies <- study_db2 %>% filter(study_sex=="mixed sex") %>% 
  arrange(desc(num_present)) %>%
  dplyr::select(-num_tot, -study_sex, -dbID, -ATC, -class, -ATC_class) %>%
  filter(num_f >= 20 & num_m >=20)
# GSE36868 - simvastatin
# GSE31312 - rituximab
# GSE100833 - ustekinumab
# GSE10846 - rituximab
# GSE13070 - insulin
# GSE47994 ... private :/ 
# GSE117239 - ustekinumab
# GSE67784 - tafamidis; has been analyzed by sex
# GSE73174 - fingolimod
# GSE9782 - bortezomib
# GSE66702 - prednisolone
# ...

large_drug_studies2 <- large_drug_studies %>% 
  dplyr::select(-num_present) %>%
  group_by(organism, data_type, study_acc) %>%
  summarize(drug=paste(name, collapse=";"), 
            num_f=unique(num_f), 
            num_m=unique(num_m)) %>%
  inner_join(exp_metadata %>% dplyr::select(study_acc, title)) %>%
  dplyr::select(organism, data_type, study_acc, drug, title, num_f, num_m) %>%
  arrange(desc(num_f+num_m))
large_drug_studies2 %>% write_csv("data/large_drug_studies.csv")


# Tafamidis-treated and untreated V30M FAP patients, asymptomatic V30M carriers,
# and healthy, age- and sex-matched controls without TTR mutations had whole blood drawn. 
my_ds <- getGEO("GSE67784")
pDat <- pData(my_ds$GSE67784_series_matrix.txt.gz)
eDat <- exprs(my_ds$GSE67784_series_matrix.txt.gz)

pDat2 <- pDat %>% 
  dplyr::select(geo_accession, title, source_name_ch1, `age:ch1`) %>%
  as_tibble() %>%
  rename(age=`age:ch1`, sample_acc=geo_accession) %>%
  mutate(sex=case_when(
    str_detect(source_name_ch1, "Female") ~ "female",
    str_detect(source_name_ch1, "Male") ~ "male")) %>%
  mutate(exposure=str_replace_all(source_name_ch1, "Male|Female| ", "")) %>%
  dplyr::select(-source_name_ch1) %>%
  mutate(sex=as.factor(sex), exposure=as.factor(exposure)) %>%
  mutate(age=as.numeric(age))
summary(pDat2)

# Symptomatic vs TafamidisTreated

pDat3 <- pDat2 %>% 
  filter(exposure %in% c("Symptomatic", "Control")) %>%
  mutate(exposure=as.factor(as.character(exposure)))
table(pDat3$sex, pDat3$exposure)

ggplot(pDat3 %>%
         mutate(grp=paste(sex,exposure)),
       aes(x=grp, y=age))+
  geom_boxplot()

# check that it matches - exposure does not
aov_test <- aov(age ~ sex + exposure + sex*exposure, data=pDat3)
summary(aov_test)

# double check sex labels
pDat4 <- pDat3 %>% 
  inner_join(comb_metadata %>%
               dplyr::select(sample_acc, metadata_sex, expr_sex, p_male),
             by=c("sample_acc")) %>%
  rename(pdat_sex=sex) %>%
  mutate(exposure=fct_recode(exposure, "control"="Symptomatic", 
                             "treated"="Control")) %>%
  mutate(pdat_sex=as.character(pdat_sex))

table(pDat4$expr_sex, pDat4$pdat_sex)
# two samples are mislabeled
pDat4 %>% filter(expr_sex!=pdat_sex) # two are labeled female but almost definitely MALE
# --> exclude these. what if they're actually other conditions!

pDat5 <- pDat4 %>% 
  filter(expr_sex==pdat_sex) %>% 
  dplyr::select(-metadata_sex, -pdat_sex, -p_male) %>%
  rename(sex=expr_sex)

#eDat2 <- eDat[,pDat5$sample_acc]
exp_dat <- read_tsv("~/Downloads/0a280a28-c929-4376-b0f8-4fd0ccb1cf15/GSE67784/GSE67784.tsv")
exp_dat2 <- data.frame(apply(exp_dat[,pDat5$sample_acc], c(1,2), function(x) log(x+1)+3))
rownames(exp_dat2) <- exp_dat$Gene

# --- attempt to use NANOVA, not worth it --- #
# f1=rep(1:2, each=8)
# f2=rep(c(1,2,1,2), each=4)
# tp=rep(1:4, 4)
# data=matrix(rnorm(16*1000), nrow=1000, ncol=16)
# result=NANOVA.test(data,f1,f2, type=1)
# fdr <-fdr.table(result)
#      c<-sig.number(fdr, FDR=0.05, qt=-1)
#result2=tanova(data,f1,f2, tp=0, longitudinal = FALSE) # fails
#result2=NANOVA.test3(data,f1,f2, tp=0, type=0) # fails
#result2=gene.classifier2(data,f1,f2, time.course=FALSE)  # fails, 1+2 fail too


res <- NANOVA.test(exp_dat2, 
              f1=as.numeric(as.factor(pDat5$sex)), 
              f2=as.numeric(as.factor(pDat5$exposure)), 
              type=1)
sigp <- which(res$pvalue < 0.01)
res$delta[sigp]
res$p[sigp]
# // TODO: add pretty names
exp_dat3<- exp_dat2[res$gene.order[1:length(sigp)],]
genes <- rownames(exp_dat3)
#rot_df <- t(exp_dat2[sigp,]) %>% as_tibble() %>% pivot_longer(cols=everything(), names_to="Gene", values_to="value")


# ---- Code for making a pretty figure ---- #
df <- cbind(t(exp_dat2[sigp,]), "sex"=pDat5$sex, "exposure"=pDat5$exposure) %>% as_tibble()
df2 <- df %>% 
  pivot_longer(ENSG00000010626:ENSG00000278685,
               names_to="gene", values_to="value")
ensembl=useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl")
queryResults <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                      filters="ensembl_gene_id",
                      values=genes,
                      mart=ensembl)
df3 <- df2 %>% left_join(queryResults, by=c("gene"="ensembl_gene_id"))
ggplot(df3 %>% mutate(grp=paste(sex, exposure)) %>%
         mutate(value=as.numeric(value)) %>%
         filter(gene %in% genes[1:8]), 
       aes(x=exposure, y=value, col=sex))+
  geom_violin(trim=FALSE, adjust=1.5) + 
  geom_point(pch="-", size=4, alpha=0.5) +
  stat_summary(fun=mean, geom="point", size=2, alpha=0.7) + 
  stat_summary(fun=mean, geom="line", size=1, 
               mapping=aes(group=sex), linetype="dashed",
               col="gray") +
  facet_grid(hgnc_symbol~sex, scales="free")+
  theme_bw()+
  theme(strip.text.y.right = element_text(angle = 0))
ggsave("figures/paper_figs/test_drug_sig_comp.png")

# ----- ------ #

fdr <-fdr.table(res)
c <- sig.number(fdr, FDR=0.05, qt=-1)

# -- try rank anova -- #
raov.full <- do.call(rbind, lapply(1:nrow(exp_dat2), function(idx){
  raov.res <- raov(gene ~ exposure*sex, data=tibble("exposure"=pDat5$exposure,
                                      "sex"=pDat5$sex,
                                     "gene"= unlist(exp_dat2[idx,])))
  raov_df <- raov.res$table %>% as_tibble()
  raov_df$term <- unlist(rownames(raov.res$table))
  raov_df$gene <- rep(rownames(exp_dat2)[idx], 3)
  return(raov_df)
}))

r2 <- raov.full %>% 
  filter(term=="exposure:sex") %>% 
  arrange(`p-value`) 
r2 %>% pull(p)
fdr.es <- fdrtool(r2 %>% pull(`p-value`), statistic="pvalue")
str(fdr.es,2)


# ---- use limma ---- #

design <- model.matrix(~exposure*sex, data=pDat5)

fit <- lmFit(exp_dat2, design)
#cont.matrix <- cbind(TvsCinF=c(0,0,-2,-2),TvsCinM=c(0,0,-2,2),Diff=c(0,0,0,4))
#fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit)

topTable(fit2)
intTop <-topTable(fit2, coef="exposurecontrol:sexmale", n=100)
head(intTop)
