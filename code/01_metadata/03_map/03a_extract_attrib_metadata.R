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


# -- useful funs -- #
common_col <- function(dat, col) {
  dat %>% 
    group_by({{col}}) %>% 
    count() %>% 
    arrange(desc(n)) %>%
    ungroup()
}
print_all <- function(dat) print(dat, n=nrow(dat))

clean_str <-function(my_str) {
  my_str %>%
    str_replace_all(PUNCT.STR , " ") %>%
    str_squish() %>%
    str_trim() %>%
    tolower()
}
####


# --- load all the data --- #
ae_attrib <- read_csv("data/ae_attributes.csv")
gsm_attrib <- read_csv("data/gsm_key_value.csv")

load("data/sample_to_attr_sm.RData") # 137 MB --> sample_dat3
sample_dat4 <- sample_dat3 %>% 
  filter(key!="biomaterial provider")

run_to_sample <- read_csv("data/sra_run_to_sample.csv")
rnaseq_dat <- sample_dat4 %>% 
  inner_join(run_to_sample, by=c("sample"="samples")) %>%
  rename(sample_acc=run) %>%
  select(sample_acc, key, value)

# --- put it all together --- #

all_attrib <- rnaseq_dat %>%
  bind_rows(gsm_attrib %>% rename(sample_acc=gsm)) %>%
  bind_rows(ae_attrib) 

# NOTE - slow
all_attrib_clean <- all_attrib %>%
  mutate(across(c(key, value), clean_str, .names="{col}_clean")) 
all_attrib_clean <- all_attrib_clean %>% select(-study_acc)



# --- 1. sex labels --- #
sg_sample_attr <- all_attrib_clean %>% 
  filter(str_detect(key_clean, "\\bsex\\b|\\bgender\\b"))
not_sg <- all_attrib_clean %>% 
  anti_join(sg_sample_attr, by="sample_acc") 
rescue_sg <- not_sg  %>% 
  filter(!str_detect(key_clean, "strain|genetic")) %>% 
  filter(str_detect(value_clean, "male"))

# sample, key, value, type_key (std/rescue), normalized_value
sg_mapped <- sg_sample_attr %>% map_sex_metadata(value)
sg_mapped2 <- rescue_sg %>% map_sex_metadata(value)
sg_tab <- sg_sample_attr %>% 
  mutate(type_key="std") %>%
  bind_rows(rescue_sg %>% mutate(type_key="rescue")) 
sg_tab2 <- sg_tab %>% 
  left_join(sg_tab %>% map_sex_metadata(value), 
                                by=c("value"="sex"))

sg_tab2 %>% write_csv("data/01_metadata/mapped_sl_all.csv")


# --- 2. cell line ---- #

# A) using keys
cell_data <- all_attrib_clean %>% 
  filter((str_detect(key_clean, "\\bcell") |
            str_detect(value_clean, "\\bcell")) & 
           !str_detect(key_clean, "tissue")) 

cell_line <- all_attrib_clean %>% 
  filter(str_detect(key_clean, "cell"),
         str_detect(key_clean, "line"))

cell_line2 <- all_attrib_clean %>% 
  anti_join(cell_line) %>%
  filter(str_detect(value_clean, "cell"),
         str_detect(value_clean, "line"))

# B) try exact match to values
map1 <- mapTextCl(cell_line %>% 
                    rename(str=value) %>% 
                    mutate(orig_str=str) %>% 
                    select(sample, orig_str, str), 
                  cell_df_nodash)
nrow(map1)
nrow(cell_line %>% distinct(value)) # 3655
map_input2 <- all_attrib_clean %>% 
  anti_join(cell_line) %>% 
  rename(str=value) %>% 
  mutate(orig_str=str) %>%
  distinct(sample, orig_str, str)
map2 <- map_input2 %>% mapTextCl(cell_df_nodash, three_l=FALSE)
test2 <- all_attrib_clean %>% 
  anti_join(cell_line) %>% 
  inner_join(map2, by=c("value"="orig_str")) 
# //TODO: age: p100, p162, etc is bad
# key: tissue, value:lncap
# 129s6 is bad -- need a list of common mouse strains to remove
# // TODO read in list of mouse strains


# --- 3. sample type --- #
all_attrib_clean2 <- all_attrib_clean %>% 
  select(-key, -value) %>% 
  rename(key=key_clean, value=value_clean)
xenografts <- all_attrib_clean2 %>% 
  filter(str_detect(key, "xenograft") | 
           str_detect(value, "xenograft"))

stem_cells <- all_attrib_clean2 %>%
  filter(str_detect(key, "\\bipsc\\b") | str_detect(value, "\\bipsc\\b") |
           str_detect(key, "\\bes cell") | str_detect(value, "\\bes cell") |
           str_detect(key, "\\besc\\b") | str_detect(value, "\\besc\\b") |
           str_detect(key, "\\bips cell") | str_detect(value, "\\bips cell") |
           ((str_detect(value, "stem") | str_detect(key, "stem")) & 
              (str_detect(key, "cell") | str_detect(value, "cell")))) 

primary_cells <- all_attrib_clean2 %>% 
  filter(str_detect(value, "primary") | str_detect(key, "primary"), 
         str_detect(key, "cell") | str_detect(value, "cell"))  

# tissue
tissue_dat <- all_attrib_clean2 %>% filter(str_detect(key, "tissue|organ") &
                                       !str_detect(key, "organism|cell|organoid")) 
rescue_tiss <- all_attrib_clean2 %>% anti_join(tissue_dat) %>%
  filter(str_detect(value, "tissue|organ") &
           !str_detect(value, "organism|cell|organoid") &
           key!="biomaterial_provider") 
# Note: "cell line: primary tissue" or "cell line: tissue"
tissue <- bind_rows(tissue_dat, rescue_tiss)


cancer_key_words <- str_squish("tumor|cancer|neoplasm|malignant|carcinoma\\b|sarcoma\\b|melanoma\\b|adenoma\\b|lymphoma\\b|leukemia\\b|mesothelioma\\b|hemangioma\\b|glioma\\b|blastoma\\b")
cancer <- all_attrib_clean2 %>% 
  filter(str_detect(key, cancer_key_words) |str_detect(value, cancer_key_words))


hc_cells <- all_attrib_clean2 %>% 
  filter(str_detect(key, "culture|passage") |
           str_detect(value, "culture|passage|\\bcell"))



cell_line_from_tiss <- all_attrib_clean2 %>% 
  filter(key=="tissue",
         str_detect(value, "\\bcell"), 
         str_detect(value, "line"))
tiss_from_cell_line <- all_attrib_clean2 %>%
  filter(str_detect(key,"\\bcell\\b"), 
         key!="progenitor cell type",
         str_detect(value, "tissue"),
         !str_detect(value, "\\bcell"))
blood <- all_attrib_clean2 %>%
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
trt_dat <- all_attrib_clean %>% 
  filter(str_detect(key, "treatment|treated|drug|compound") |
           str_detect(value, "treatment|treated|drug|compound")) 
trt_in <- trt_dat %>% distinct(value) %>% rename(str=value) %>% mutate(src_col=str)
drug_dat <- labelNgram(trt_in, drug_info_df)

# B) exact match to values

trt_dat2 <- all_attrib_clean %>% 
  anti_join(trt_dat)

trt_in2 <- trt_dat2 %>% 
  distinct(value) %>% 
  rename(str=value) %>% 
  mutate(src_col=str)
drug_dat2 <- labelNgram(trt_in2, drug_info_df)

# same, "olive oil", "nadh", "dmso", "lactose", "dhea", "cyclo", "beam",
