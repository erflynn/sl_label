# 
# TODO:
# - reorganize files!!!
# - add in mouse sample strain
# - add in list of major tissue types
# - get arrayexpress data
# - check GEO data
# - clean up CL (needs Ngram) + drug mapping functions
# - apply cleaned functions --> drugs, cells
# - create output df for sample source
# - make sure using "\\b" when we need

require('tidyverse')
source("code/01_metadata/03_map/00_mapping_utils.R")
source("code/01_metadata/03_map/00_map_sex_metadata_f.R")

ae_attrib <- read_csv("data/ae_attributes.csv")
gsm_attrib <- read_csv("data/gsm_key_value.csv")

# -- useful funs -- #
common_col <- function(dat, col) {
  dat %>% 
    group_by({{col}}) %>% 
    count() %>% 
    arrange(desc(n)) %>%
    ungroup()
}
print_all <- function(dat) print(dat, n=nrow(dat))
####

load("data/sample_to_attr_sm.RData") # 137 MB --> sample_dat3
sample_dat4 <- sample_dat3 %>% 
  mutate(key=str_trim(str_squish(str_replace_all(key, PUNCT.STR , " ")))) %>%
  mutate(value=str_trim(str_squish(str_replace_all(value, PUNCT.STR , " ")))) %>%
  filter(key!="biomaterial provider")



# --- 1. sex labels --- #
sg_sample_attr <- sample_dat3 %>% 
  filter(str_detect(key, "sex|gender") & ! 
           str_detect(key, "sexually transmitted"))
not_sg <- sample_dat3 %>% anti_join(sg_sample_attr, by="sample") 
rescue_sg <- not_sg  %>% 
  filter(!str_detect(key, "strain|genetic")) %>% 
  filter(str_detect(value, "male"))

# sample, key, value, type_key (std/rescue), normalized_value
sg_mapped <- sg_sample_attr %>% map_sex_metadata(value)
sg_mapped2 <- rescue_sg %>% map_sex_metadata(value)
sg_tab <- sg_sample_attr %>% 
  mutate(type_key="std") %>%
  bind_rows(rescue_sg %>% mutate(type_key="rescue")) 
sg_tab2 <- sg_tab %>% 
  left_join(sg_tab %>% map_sex_metadata(value), 
                                by=c("value"="sex"))


sg_tab2 %>% write_csv("data/01_metadata/mapped_rnaseq_sex.csv")

# collapse
sg_tab3 <- sg_tab2 %>% 
  filter(mapped_sex!="") %>%
  group_by(sample) %>%
  summarise(mapped_sex=paste(unique(mapped_sex), collapse=";"),
            type_key=paste(unique(type_key), collapse=";")) 
sg_tab4 <- sg_tab3 %>%
  mutate(mapped_sex=ifelse(str_detect(mapped_sex, ";"),"unknown", mapped_sex))
stopifnot((sg_tab4 %>% distinct(sample) %>% nrow())==nrow(sg_tab4))            

# add a run column
sg_run <- run_to_sample %>% inner_join(sg_tab4, by=c("samples"="sample")) %>%
  rename(sample=samples, sex=mapped_sex, key_type=type_key) 
stopifnot((sg_run %>% distinct(run) %>% nrow())==nrow(sg_run))            
sg_run %>% write_csv("data/01_metadata/run_mapped_rnaseq_sex.csv")

# compare
comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", 
                          col_types="cccccccdcd")

rnaseq_dat <- comb_metadata %>% filter(data_type=="rnaseq")
rnaseq2 <- rnaseq_dat %>% 
  select(sample_acc, organism, metadata_sex) %>% 
  left_join(sg_run, by=c("sample_acc"="run")) %>%
  mutate(sex=ifelse(is.na(sex), "unknown", sex))
table(rnaseq2$metadata_sex, rnaseq2$sex)

sum(rnaseq2$metadata_sex=="unknown")-sum(rnaseq2$sex=="unknown") # 45.7k
# better :)

# not organism-specific
rnaseq2 %>% 
  group_by(organism) %>% 
  summarise(ms=sum(metadata_sex=="unknown"), ss=sum(sex=="unknown"), n=n()) %>%
  mutate(diff=ms-ss) %>%
  mutate(across(c(ss, ms, diff), ~./n))

#sum(rnaseq2$sex=="unknown")/nrow(rnaseq2) # 0.7600228
#sum(rnaseq2$metadata_sex=="unknown")/nrow(rnaseq2) # 0.8376703

rnaseq2 %>% 
  filter(metadata_sex=="unknown" & sex!="unknown") %>%
  common_col(key_type) # 44482 come from std, 1373 come from rescue

rnaseq2 %>% 
  filter(metadata_sex=="unknown" & sex!="unknown") %>%
  sample_n(10) %>%
  left_join(sg_tab2, by="sample") %>% select(key, value, sex)

# compare to metasra/phepred
load("data/10_comparison/metasra_data.RData")

metasra_sl <- mapped_df %>% 
  filter(term_id %in% c("UBERON:0003100", "UBERON:0003101")) %>%
  mutate(metasra_sex=case_when(term_id=="UBERON:0003100" ~ "female",
                               term_id=="UBERON:0003101" ~ "male"))  %>%
  group_by(sample_accession) %>%
  summarise(metasra_sex=paste(unique(metasra_sex), collapse=";")) %>%
  mutate(metasra_sex=ifelse(str_detect(metasra_sex, ";"), "mixed", metasra_sex)) %>%
  ungroup() 
nrow(metasra_sl) # 153417
nrow(sg_tab4) # 127755

# NOTE - metasra has no mouse! but somehow has more data? but less for our samples

comparison_tab <- metasra_sl %>%   
  left_join(run_to_sample, by=c("sample_accession"="samples")) %>%
  right_join(rnaseq2, by=c("run"="sample_acc")) %>%
  mutate(metasra_sex=ifelse(is.na(metasra_sex), "unknown", metasra_sex)) %>%
  rename(rb_sex=metadata_sex, mapped_sex=sex, sample_acc=run) %>%
  select(-sample_accession, -sample,-key_type) %>%
  select(sample_acc, organism, everything())
comparison_tab %>% 
  filter(organism=="human") %>% 
  summarise(ms=sum(metasra_sex=="unknown"), ss=sum(mapped_sex=="unknown"), n=n()) %>%
  mutate(diff=ms-ss) %>%
  mutate(across(c(ss, ms, diff), ~./n))

conf_mat_human <- comparison_tab %>% 
  filter(organism=="human") %>% 
  group_by(metasra_sex, mapped_sex) %>% 
  count() %>% 
  ungroup() %>%
  pivot_wider(names_from="mapped_sex", values_from="n", values_fill=0) %>%
  select(metasra_sex, female, male, mixed, unknown)
# columns are our labels
conf_mat_human %>% write_csv("tables/metasra_conf_mat_sex.csv")

# --- 2. cell line ---- #

# A) using keys
cell_data <- sample_dat4 %>% 
  filter((str_detect(key, "\\bcell") | str_detect(value, "\\bcell")) & 
           !str_detect(key, "tissue")) 


cell_line <- sample_dat4 %>% 
  filter(str_detect(key, "cell") & str_detect(key, "line"))

cell_line2 <- sample_dat4 %>% anti_join(cell_line) %>%
  filter(str_detect(value, "cell") & str_detect(value, "line"))

# B) try exact match to values
map1 <- mapTextCl(cell_line %>% 
                    rename(str=value) %>% 
                    mutate(orig_str=str) %>% 
                    select(sample, orig_str, str), 
                  cell_df_nodash)
nrow(map1)
nrow(cell_line %>% distinct(value)) # 3655
map_input2 <- sample_dat4 %>% 
  anti_join(cell_line) %>% 
  rename(str=value) %>% 
  mutate(orig_str=str) %>%
  distinct(sample, orig_str, str)
map2 <- map_input2 %>% mapTextCl(cell_df_nodash, three_l=FALSE)
test2 <- sample_dat4 %>% 
  anti_join(cell_line) %>% 
  inner_join(map2, by=c("value"="orig_str")) 
# //TODO: age: p100, p162, etc is bad
# key: tissue, value:lncap
# 129s6 is bad -- need a list of common mouse strains to remove
# // TODO read in list of mouse strains


# --- 3. sample type --- #
xenografts <- sample_dat4 %>% 
  filter(str_detect(key, "xenograft") | str_detect(value, "xenograft"))

stem_cells <- sample_dat4 %>%
  filter(str_detect(key, "\\bipsc\\b") | str_detect(value, "\\bipsc\\b") |
           str_detect(key, "\\bes cell") | str_detect(value, "\\bes cell") |
           str_detect(key, "\\besc\\b") | str_detect(value, "\\besc\\b") |
           str_detect(key, "\\bips cell") | str_detect(value, "\\bips cell") |
           ((str_detect(value, "stem") | str_detect(key, "stem")) & 
              (str_detect(key, "cell") | str_detect(value, "cell")))) 

primary_cells <- sample_dat4 %>% 
  filter(str_detect(value, "primary") | str_detect(key, "primary"), 
         str_detect(key, "cell") | str_detect(value, "cell"))  

# tissue
tissue_dat <- sample_dat4 %>% filter(str_detect(key, "tissue|organ") &
                                       !str_detect(key, "organism|cell|organoid")) 
rescue_tiss <- sample_dat4 %>% anti_join(tissue_dat) %>%
  filter(str_detect(value, "tissue|organ") &
           !str_detect(value, "organism|cell|organoid") &
           key!="biomaterial_provider") 
# Note: "cell line: primary tissue" or "cell line: tissue"
tissue <- bind_rows(tissue_dat, rescue_tiss)


cancer_key_words <- str_squish("tumor|cancer|neoplasm|malignant|carcinoma\\b|sarcoma\\b|melanoma\\b|adenoma\\b|lymphoma\\b|leukemia\\b|mesothelioma\\b|hemangioma\\b|glioma\\b|blastoma\\b")
cancer <- sample_dat4 %>% 
  filter(str_detect(key, cancer_key_words) |str_detect(value, cancer_key_words))


hc_cells <- sample_dat4 %>% 
  filter(str_detect(key, "culture|passage") |
           str_detect(value, "culture|passage|\\bcell"))



cell_line_from_tiss <- sample_dat4 %>% 
  filter(key=="tissue",
         str_detect(value, "\\bcell"), 
         str_detect(value, "line"))
tiss_from_cell_line <- sample_dat4 %>%
  filter(str_detect(key,"\\bcell\\b"), 
         key!="progenitor cell type",
         str_detect(value, "tissue"),
         !str_detect(value, "\\bcell"))
blood <- sample_dat4 %>%
  filter(str_detect(key, "pbmc|whole blood") | 
           str_detect(value,"pbmc|whole blood" ))

# sample_id | xenograft | stem_cells | primary_cells | tissue | cancer | cell_line | cells
sample_src_tab <- xenografts %>% select(sample) %>% mutate(xenograft=TRUE) %>%
  full_join(stem_cells %>% select(sample) %>% mutate(stem_cell=TRUE)) %>%
  full_join(primary_cells %>% select(sample) %>% mutate(primary_cell=TRUE)) %>%
  full_join(tissue %>% select(sample) %>% mutate(tissue=TRUE)) %>%
  full_join(cancer %>% select(sample) %>% mutate(cancer=TRUE)) %>%
  full_join(cell_line %>% bind_rows(cell_line2) %>% select(sample) %>% mutate(cell_line=TRUE)) %>%
  full_join(hc_cells %>% select(sample) %>% mutate(cells=TRUE))

sample_src_tab2 <- sample_src_tab %>%
  mutate(across(-sample, ~replace_na(.,FALSE))) %>%
  distinct()
stopifnot(length(unique(sample_src_tab2$sample))==nrow(sample_src_tab2))

sample_src_tab2 %>% filter(primary_cell, tissue)

src_tab <- sample_src_tab2 %>%
  mutate(source_type=case_when(
    xenograft ~ "xenograft",
    stem_cell ~ "stem cells",
    primary_cell~ "primary cells",
    cell_line ~ "cell line",
    cells ~ "cell",
    tissue ~ "tissue",
    #cells ~ "cell",
    TRUE~"other"
  )) %>%
  mutate(source_type=ifelse(sample %in% blood$sample | 
                               sample %in% tiss_from_cell_line$sample, 
                            "tissue", source_type))

table(src_tab$source_type)

src_tab %>% filter(tissue, cell_line) %>% nrow()

compare_src <- src_tab %>% 
  select(sample, source_type) %>%
  full_join(run_to_sample, by=c("sample"="samples")) %>%
  mutate(source_type=replace_na(source_type, "other")) %>%
  left_join(sample_type_df %>% rename(metasra_type=sample_type), 
            by=c("sample"="sample_accession")) %>%
  mutate(metasra_type=case_when(is.na(metasra_type) ~ "other",
                                str_detect(metasra_type, "induced") ~ "stem cells",
                                TRUE ~ metasra_type) )
  
compare_src %>% 
  mutate(metasra_type=factor(metasra_type, levels=c("in vitro differentiated cells","cell line", "primary cells", 
                                                    "stem cells", "tissue", "other"
                                                    )),
         source_type=factor(source_type, levels=c("cell", "cell line", "primary cells", 
                                                   "stem cells", "tissue", "other",
                                                    "xenograft"))) %>%
  group_by(metasra_type, source_type) %>%
  count() %>%
  pivot_wider(names_from="source_type", values_from="n", values_fill=0)
compare_src %>% 
  filter(!metasra_type %in% c( "other" , "in vitro differentiated cells") &
           !source_type %in% c("other", "xenograft", "cell")) %>%
  mutate(match=(metasra_type==source_type)) %>%
  group_by(match) %>% count() 

simple_comp <- compare_src %>% 
  mutate(
    metasra_type=case_when(
      metasra_type=="other" ~ "other",
      metasra_type=="tissue" ~ "tissue",
      str_detect(metasra_type, "cell")  ~ "cell"),
    source_type=case_when(
      source_type=="other" ~ "other",
      source_type=="tissue" ~ "tissue",
      str_detect(source_type, "cell")  ~ "cell",
      source_type=="xenograft" ~ "xenograft"
    )) %>%
  ungroup()
simple_comp %>%
  group_by(metasra_type, source_type) %>% count() %>%
  pivot_wider(names_from=source_type, values_from=n) 
simple_comp %>% 
  filter(metasra_type %in% c("cell", "tissue", "other"),
         source_type %in% c("cell", "tissue", "other")) %>%
  mutate(across(c(metasra_type, source_type), 
                ~factor(.,levels=c("cell", "tissue", "other")))) %>%
  mutate(tot=n()) %>%
  group_by(metasra_type, source_type, tot) %>% count() %>%
  pivot_wider(names_from=source_type, values_from=n) %>%
  mutate(across(c(cell, tissue, other), ~./tot))

simple_comp %>% 
  filter(metasra_type=="cell" & source_type=="tissue") %>%
  filter(confidence > 0.9) %>%
  sample_n(10)

simple_comp %>% 
  filter(metasra_type=="tissue" & source_type=="other") %>%
  filter(confidence > 0.9) %>%
  sample_n(10)

# TODO compare to JSON http://metasra.biostat.wisc.edu/publication.html

# ---- 4. Drugs--- #
# A) keys
trt_dat <- sample_dat4 %>% 
  filter(str_detect(key, "treatment|treated|drug|compound") |
           str_detect(value, "treatment|treated|drug|compound")) 
trt_in <- trt_dat %>% distinct(value) %>% rename(str=value) %>% mutate(src_col=str)
drug_dat <- labelNgram(trt_in, drug_info_df)

# B) exact match to values

trt_dat2 <- sample_dat4 %>% 
  anti_join(trt_dat)

trt_in2 <- trt_dat2 %>% 
  distinct(value) %>% 
  rename(str=value) %>% 
  mutate(src_col=str)
drug_dat2 <- labelNgram(trt_in2, drug_info_df)

# same, "olive oil", "nadh", "dmso", "lactose", "dhea", "cyclo", "beam",
