# Look at sex differences in drug response in studies that are sufficiently powered
#
#

require('GEOquery')
require('limma')
require('tidyverse')

# ---- are there any studies where we could look at sex diff in drug response? ---- $
study_db2 %>% filter(study_sex=="mixed sex") %>% arrange(desc(num_present)) %>%
  select(-num_tot, -study_sex, -dbID, -ATC, -class, -ATC_class) %>%
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

# Tafamidis-treated and untreated V30M FAP patients, asymptomatic V30M carriers,
# and healthy, age- and sex-matched controls without TTR mutations had whole blood drawn. 
my_ds <- getGEO("GSE67784")
pDat <- pData(my_ds$GSE67784_series_matrix.txt.gz)
eDat <- exprs(my_ds$GSE67784_series_matrix.txt.gz)

pDat2 <- pDat %>% 
  select(geo_accession, title, source_name_ch1, `age:ch1`) %>%
  as_tibble() %>%
  rename(age=`age:ch1`, sample_acc=geo_accession) %>%
  mutate(sex=case_when(
    str_detect(source_name_ch1, "Female") ~ "female",
    str_detect(source_name_ch1, "Male") ~ "male")) %>%
  mutate(exposure=str_replace_all(source_name_ch1, "Male|Female| ", "")) %>%
  select(-source_name_ch1) %>%
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
               select(sample_acc, metadata_sex, expr_sex, p_male),
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
  select(-metadata_sex, -pdat_sex, -p_male) %>%
  rename(sex=expr_sex)

#eDat2 <- eDat[,pDat5$sample_acc]
exp_dat <- read_tsv("~/Downloads/0a280a28-c929-4376-b0f8-4fd0ccb1cf15/GSE67784/GSE67784.tsv")
exp_dat2 <- data.frame(apply(exp_dat[,pDat5$sample_acc], c(1,2), function(x) log(x+1)+3))
rownames(exp_dat2) <- exp_dat$Gene




design <- model.matrix(~exposure*sex, data=pDat5)
fit <- lmFit(exp_dat2, design)
#cont.matrix <- cbind(TvsCinF=c(0,0,-2,-2),TvsCinM=c(0,0,-2,2),Diff=c(0,0,0,4))
#fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit)

topTable(fit2)
topTable(fit2, coef="exposurecontrol:sexmale", )
# -- gene ~ exposure + age + sex + exposure*sex -- #
# gene ~ exp + age + exp*age + exp*sex #

