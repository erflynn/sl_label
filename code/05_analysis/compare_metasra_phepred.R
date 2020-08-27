# download the MetaSRA and phenopredict mappings!

# how does it work to map these to them?
# what does the breakdown look lik here

require('RSQLite')
require('SRAdb')
require('recount') ## warning there are multiple recount packages! --> had to remove.packages() for the original


##  ---- phenopredict ---- ##
phepred <- recount::add_predictions() # this downloads everything. 70,479
save(phepred, file="data/10_comparison/phepred_data.RData")
phepred_sl <- phepred %>% 
  select(sample_id, reported_sex, predicted_sex) %>%
  mutate(reported_sex=tolower(as.character(reported_sex)))

phepred_sl2 <- phepred_sl %>%
  mutate(reported_mapped=case_when(
    is.na(reported_sex) ~ "unknown",
    reported_sex %in% c("mixed", "mixture") ~ "mixed",
    reported_sex %in% c("male", "m") ~ "male",
    reported_sex %in% c("female", "f") ~ "female",
    reported_sex %in% c("unknown", "dk", "u", "not determined", "n/a", 
                        "not collected", "not applicable", "missing", "not available",
                        "asexual") ~ "unknown",
    str_detect(reported_sex, "male and female") ~ "mixed",
    str_detect(reported_sex, "female") & str_detect(reported_sex, " male") ~ "mixed",
    str_detect(reported_sex, "female") ~ "female",
    TRUE ~ "unknown"
  ))

## ---- metasra ---- #
con <- dbConnect(RSQLite::SQLite(), "../metasra.v1-6.sqlite")
#dbListTables(con)
sample_type_df <- dbGetQuery(con, "SELECT * FROM sample_type;") # 338,842
dbDisconnect(con)

mapped_df <- dbGetQuery(con, "SELECT * FROM mapped_ontology_terms;") # 1,858,943
rv_df <- dbGetQuery(con, "SELECT * FROM real_value_properties;") # 125,851

# map to sample IDs
sracon <- dbConnect(RSQLite::SQLite(),"../SRAmetadb.sqlite")
samp_to_runs <- sraConvert(sample_type_df$sample_accession, 'run', sracon) # 476,416  
#length(unique(samp_to_runs$run)) # 443104 --> some runs map to more than one sample!! :/
dbDisconnect(sracon)

save(mapped_df, rv_df, sample_type_df, samp_to_runs, file="data/10_comparison/metasra_data.RData")

mapped_df2 <- mapped_df %>% 
  left_join(samp_to_runs, by=c("sample_accession"="sample")) %>%
  dplyr::rename(acc=run)


unique(rv_df$property_term_id)
# "EFO:0000246" - age
# "EFO:0007061" - passage number
# "EFO:0000724" - timepoint
# "EFO:0000721" - time
# "EFO:0004340" - BMI
# "EFO:0005056" - age at death
# "EFO:0004918" - age at diagnosis
# "EFO:0000410" - disease staging

#mapped_df %>% filter(sample_accession=="SRS5483655")
#"UBERON:0003101" <- male organism "UBERON:0003100" <- "female organism"
#"UBERON:0002107" <- liver
# rv_df %>% filter(sample_accession=="SRS5483655") 
# "EFO:0000246" - age 34
# sample_type_df %>% filter(sample_accession=="SRS5483655")


# ---- clean up the labels ---- #
metasra_sl <- mapped_df %>% 
  filter(term_id %in% c("UBERON:0003100", "UBERON:0003101")) %>%
  mutate(metasra_sex=case_when(term_id=="UBERON:0003100" ~ "female",
                               term_id=="UBERON:0003101" ~ "male")) %>%
  left_join(samp_to_runs, by=c("sample_accession"="sample")) %>% 
  rename(acc=run)


# ---- what is the overlap with refine-bio data? ---- #
rnaseq_mapping <- read_csv("data/01_metadata/mouse_rnaseq_exp_to_sample.csv") # 305642
#### THERE IS SOMETHING WRONG. :(
count_dat <- read_csv("data/01_metadata/human_exp_to_sample_counts.csv")
count_dat2 <- count_dat %>% filter(present & num_reads > 100000) # 117,032

metasra_metadata_overlap <- length(intersect(count_dat$sample_acc, mapped_df2$acc)) # 66,961
metasra_data_overlap <- length(intersect(count_dat2$sample_acc, mapped_df2$acc)) # 65,387

phe_metadata_overlap <- length(intersect(count_dat$sample_acc, phepred$sample_id)) # 26,928
phe_data_overlap <- length(intersect(count_dat2$sample_acc, phepred$sample_id)) # 26,264

phe_metasra_overlap <- length(intersect(phepred$sample_id, mapped_df2$acc))  # 53,370

# ---- load refine-bio labels ---- #
rb_sex <- read.csv("data/01_metadata/human_rnaseq_metadata_sex.csv", stringsAsFactors = FALSE) %>% unique()

# ---- overlapping metadata sex ---- #


# ---- overlapping imputed sex ---- #

# cell line

# output


# --------- sex labels!!! ---------- #
# Q1 : does this match refine-bio?
##
table(sex_lab$metasra_sex)
#female   male 
#84903  70228 
sex_lab2 <- sex_lab %>% 
table(sex_lab2$metasra_sex)
# female   male 
# 162124 135183 
dim(rb_sex) # 229,789

length(setdiff(samp_to_runs$run, rb_sex$acc)) # 310,396
length(setdiff(rb_sex$acc,samp_to_runs$run)) # 97,081
length(intersect(samp_to_runs$run, rb_sex$acc)) # 132,708

length(setdiff(samp_to_runs$run, count_dat$sample_acc)) # 373,245
length(setdiff(count_dat$sample_acc,samp_to_runs$run)) # 50,141


# unique(sapply((sex_lab$sample_accession), function(x) substr(x, 1,3)))

length(intersect(rb_sex$acc, sex_lab2$acc)) # 38,441
length(setdiff(sex_lab2$acc, rb_sex$acc)) # 219,461
length(setdiff(rb_sex$acc, sex_lab2$acc)) # 191,348


comb_lab <- rb_sex %>% select(acc, mapped_sex) %>% rename(rb_sex=mapped_sex) %>% 
  inner_join(sex_lab2) %>%
  mutate(metasra_sex=ifelse(is.na(metasra_sex), "unknown", metasra_sex)) %>%
  select(-term_id) %>%
  select(-sample_accession) %>%
  unique()

#length(unique(comb_lab$acc))
#[1] 38441
# nrow(comb_lab)
# 39063

# very few disagree! 
table(comb_lab[,c("rb_sex", "metasra_sex")])
# metasra_sex
# rb_sex    female  male
# female    6413     0
# male         5  5674
# mixed      184   184
# unknown  13774 12829

# ---- grab phenopredict sex labels ----- #


# Q2: do phenopredict sex labels match ours?
length(intersect(phepred_sl$sample_id, rb_sex$acc)) # 41,697
length(setdiff(phepred_sl$sample_id, rb_sex$acc)) # 28,782
length(setdiff(rb_sex$acc, phepred_sl$sample_id)) # 188,092

# load our labels
expr_lab <- read_csv("data/09_model_summary/human_rnaseq_sex_labels.csv") %>%
  mutate(exprsex=round(pred)) 
length(intersect(rb_sex$acc, expr_lab$id)) # 116,716
length(setdiff(rb_sex$acc, expr_lab$id)) # 113,073 (expr_lab is a subset)
length(setdiff(expr_lab$id,rb_sex$acc)) # 0

my_sl <- rb_sex %>% 
  select(acc, mapped_sex) %>% 
  dplyr::rename(metadata=mapped_sex) %>% 
  left_join(expr_lab, by=c("acc"="id")) %>%
  mutate(exprsex=case_when(is.na(pred) ~ "unknown",
                           exprsex==0 ~ "female", 
                           exprsex==1 ~ "male"))



phepred_comp <- inner_join(phepred_sl2, my_sl, by=c("sample_id"="acc")) 
table(phepred_comp[,c("reported_mapped", "metadata")])
table(phepred_comp[,c("predicted_sex", "exprsex")]) 
# exprsex
#predicted_sex female  male unknown
#female      13935   691    7102
#male          755  9969    6182
#Unassigned    760    66    2237
phepred_comp %>% 
  filter(predicted_sex %in% c("female", "male") & 
           exprsex %in% c("female", "male")) %>%
  group_by(exprsex, predicted_sex) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(eq=(exprsex==predicted_sex)) %>%
  group_by(eq) %>% 
  summarize(n=sum(n)) %>% 
  ungroup() %>%
  mutate(tot=sum(n)) %>% 
  group_by(eq) %>%
  mutate(frac=n/tot) # 94.3% concordance...

# no labels -- still looks good!
phepred_comp %>% 
  filter(predicted_sex %in% c("female", "male") & 
           exprsex %in% c("female", "male") &
           reported_mapped=="unknown" &
           metadata=="unknown") %>%
  group_by(exprsex, predicted_sex) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(eq=(exprsex==predicted_sex)) %>%
  group_by(eq) %>% 
  summarize(n=sum(n)) %>% 
  ungroup() %>%
  mutate(tot=sum(n)) %>% 
  group_by(eq) %>%
  mutate(frac=n/tot) # 94.4%

# which performs better? :P
consens <- phepred_comp %>% 
  mutate(consensus_meta=case_when(
    reported_mapped=="unknown" & metadata %in% c("female", "male") ~ metadata,
    metadata=="unknown" & reported_mapped %in% c("female", "male") ~ reported_mapped,
    TRUE ~ "unknown"
    ))

comp <- consens %>%  
  filter(predicted_sex %in% c("female", "male") & 
        exprsex %in% c("female", "male")) %>%
  filter(consensus_meta != "unknown") %>%
  select(sample_id, consensus_meta, predicted_sex, exprsex, pred) 


comp %>%
  group_by(exprsex, consensus_meta) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(eq=(exprsex==consensus_meta)) %>%
  group_by(eq) %>% 
  summarize(n=sum(n)) %>% 
  ungroup() %>%
  mutate(tot=sum(n)) %>% 
  group_by(eq) %>%
  mutate(frac=n/tot) # 84.6%

comp %>% 
  group_by(predicted_sex, consensus_meta) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(eq=(predicted_sex==consensus_meta)) %>%
  group_by(eq) %>% 
  summarize(n=sum(n)) %>% 
  ungroup() %>%
  mutate(tot=sum(n)) %>% 
  group_by(eq) %>%
  mutate(frac=n/tot) # 87.0%

# soooo.... we perform worse:
#   possibilities:
#     - their training data is in there
gtex <- predictions %>% filter(dataset=="gtex")
comp %>% semi_join(gtex) # gtex comparison - ok so it's not in there

# // TODO: predict on GTEX



# which do we differ on?
# -- are these cell line labels? we appear to be overpredcting f
table(predictions$predicted_samplesource)
cl_line <- predictions %>% filter(predicted_samplesource=="cell_line")
comp %>% anti_join(cl_line) %>%
  group_by(exprsex, consensus_meta) %>%
  dplyr::count() # 85.6%

comp %>% anti_join(cl_line) %>%
  group_by(predicted_sex, consensus_meta) %>%
  dplyr::count() # 88.1%

# they're still doing better

# well theyre mostly in the fuzzy range...
comp %>% filter(consensus_meta==predicted_sex & exprsex != predicted_sex)
comp %>% filter(consensus_meta==exprsex & exprsex != predicted_sex)

# what do they put as unassigned?


# --------- cell line labels ---------- #

sample_src <- read_csv("data/sample_source_type.csv") %>% select(acc,source_type)
samp_to_runs <- read_csv("data/sra_run_to_sample.csv") # 588295
run_type <- sample_type_df %>%  # 338842
  inner_join(samp_to_runs, by=c("sample_accession"="samples")) # 139432
compare_df <- run_type %>% inner_join(sample_src, by=c("run"="acc"))

compare_df2 <- compare_df %>% 
  mutate(source_type=case_when(
    source_type %in% c("unnamed_cl", "named_cl") ~ "cell line",
    source_type=="primary_cells" ~ "primary cells",
    source_type=="stem_cell" ~ "stem cells",
    source_type=="tissue" ~ "tissue",
    TRUE ~ "other")) %>%
  mutate(sample_type=ifelse(sample_type=="induced pluripotent stem cell line", "stem cells", sample_type))
table(compare_df2$sample_type, compare_df2$source_type)

# we are labeling a lot of cell line data as tissue
# we are labeling a lot of primary cells as tissue
compare_df2 %>% filter(sample_type=="cell line" & source_type=="tissue") %>%
  sample_n(50)
# our error: SRS859083

# Q3: does our group label match theirs for both MetaSRA and phenopredict (e.g. cell line, non-cell line)?
samp2 <- sample_type_df %>% 
  left_join(samp_to_runs, by=c("sample_accession"="sample")) %>% 
  dplyr::rename(acc=run) %>%
  select(-sample_accession) %>% 
  unique()
phepred_cl <- predictions %>% select(sample_id, predicted_samplesource)
length(intersect(phepred_cl$sample_id, samp2$acc)) # 56188
length(setdiff(phepred_cl$sample_id, samp2$acc)) # 14291
length(setdiff( samp2$acc, phepred_cl$sample_id)) # 386917

samp3 <- samp2 %>% mutate(sample_type2=case_when(sample_type=="cell line" ~ "cell_line",
                                        sample_type=="tissue" ~ "tissue",
                                        TRUE ~ "other")) 

comp_pred_cl2 <- phepred_cl %>% inner_join(samp3, by=c("sample_id"="acc"))
table(comp_pred_cl2[,c("predicted_samplesource", "sample_type")])
  
table(comp_pred_cl2[,c("predicted_samplesource", "sample_type2")])
comp_pred_cl2 %>% filter(sample_type2 %in% c("cell_line", "tissue") &
                           predicted_samplesource %in% c("cell_line", "tissue")) %>%
  group_by(predicted_samplesource, sample_type2) %>%
  dplyr::count() %>%
  ungroup() %>%
  mutate(eq=(predicted_samplesource==sample_type2)) %>%
  group_by(eq) %>% 
  summarize(n=sum(n)) %>% 
  ungroup() %>%
  mutate(tot=sum(n)) %>% 
  group_by(eq) %>%
  mutate(frac=n/tot) # 79.7% concordance

# what about our labels?

sample_cl <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
sample_cl2 <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl_part2.csv")
sample_cl3 <- sample_cl %>% select(gsm, accession) %>% bind_rows(sample_cl2 %>% select(gsm, accession))

# add IN/not IN labels
mapping2 <- mapping %>% anti_join(sample_cl3, by=c("sample_acc"="gsm"))
table((samp3 %>% inner_join(mapping, by=c("acc"="sample_acc")))$sample_type2)
cl_both <- samp3 %>% inner_join(sample_cl3, by=c("acc"="gsm"))
table(cl_both$sample_type2)

# Q4: does our mapped label match theirs for MetaSRA?
require('rjson')
cvcl_map <- fromJSON(file="../metasra_cvcl_mappings.json") 

cvcl_df <- do.call(rbind, lapply(cvcl_map, function(x) paste(sort(x$mapped_terms), collapse=";")))
cvcl_df %>% head()
cell_mapped <- mapped_df %>% filter(str_detect(term_id, "CVCL")) %>%
  mutate(term_id2=str_replace_all(tolower(term_id), ":", "_"))
cl_mapped2 <- cell_mapped %>% left_join(samp_to_runs, by=c("sample_accession"="sample")) %>%
  dplyr::rename(acc=run) 
both_cl <- cl_mapped2 %>% inner_join(sample_cl3, by=c("acc"="gsm"))
both_cl %>% filter(term_id2==accession) %>% nrow() # 11738/16951 (69.2%)
multi_map <- both_cl %>% filter(term_id2!=accession & str_detect(accession, ";")) 
some_overlap <- multi_map %>% filter(str_detect(accession,term_id2)) # 1962 --> 72.3% show some overlap
no_overlap <- multi_map %>% filter(!str_detect(accession,term_id2)) # 548

# // TODO - bug, repeated CL, e.g. "cvcl_iq55;cvcl_6502;cvcl_6502

both_cl %>% filter(term_id2!=accession & !str_detect(accession, ";")) # 2703 do not match
both_cl %>% filter(term_id2!=accession & !str_detect(accession, ";")) %>% select(term_id2, accession) %>% unique() %>% head()
metadata <- read.csv("data/01_metadata/human_rnaseq_sample_metadata.csv") 
metadata %>% filter(acc=="ERR2929108") # ours looks better lol

# what is the overlap??



# ------- OTHER ------- #

# sample type



