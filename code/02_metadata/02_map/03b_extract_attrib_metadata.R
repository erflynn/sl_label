# 
# TODO:
# - reorganize files!!!
# - add in mouse sample strain
# - add in list of major tissue types
# - get arrayexpress data
# - check GEO data
# - apply cleaned functions --> drugs, cells
# - create output df for sample source
# - make sure using "\\b" when we need

require('tidyverse')
source("code/01_metadata/03_map/00_mapping_utils.R")
source("code/01_metadata/03_map/00_map_sex_metadata_f.R")


# -- useful funs -- #

clean_str <-function(my_str, fill=" ") {
  my_str %>%
    str_replace_all(PUNCT.STR , fill) %>%
    str_squish() %>%
    str_trim() %>%
    tolower()
}




####
load("data/01_metadata/all_sample_attrib_clean.RData") # --> all_attrib_clean


# --- 1. sex labels --- #
sg_sample_attr <- all_attrib_clean %>% 
  filter(str_detect(key_clean, "\\bsex\\b|\\bgender\\b"))
not_sg <- all_attrib_clean %>% 
  anti_join(sg_sample_attr, by="sample_acc") 
rescue_sg <- not_sg  %>% 
  filter(!str_detect(key_clean, "strain|genetic|recipient|donor|child|sex chromosome complement")) %>% 
  filter(str_detect(value_clean, "male"))

# sample, key, value, type_key (std/rescue), normalized_value
sg_mapped <- sg_sample_attr %>% map_sex_metadata(value_clean)
sg_mapped2 <- rescue_sg %>% map_sex_metadata(value_clean)
sg_tab <- sg_sample_attr %>% 
  mutate(type_key="std") %>%
  bind_rows(rescue_sg %>% mutate(type_key="rescue")) 
sg_tab2 <- sg_tab %>% 
  left_join(sg_tab %>% map_sex_metadata(value_clean), 
                                by=c("value_clean"="sex"))


sg_tab3 <- sg_tab2 %>% unique() %>%   filter(!str_detect(key_clean, "recipient|donor|child|sex chromosome complement")) 
  
mult_row <- sg_tab3 %>%  group_by(sample_acc) %>% count() %>% arrange(desc(n)) %>% filter(n > 1)
single_row <- sg_tab3 %>% anti_join(mult_row) %>% select(sample_acc, key_clean, value_clean, type_key, mapped_sex)
mult_row_condensed <- sg_tab3 %>% semi_join(mult_row) %>%
  filter(mapped_sex != "unknown") %>%
  group_by(sample_acc) %>%
  summarise(key_clean=paste(unique(key_clean), collapse=";"),
            value_clean=paste(unique(value_clean), collapse=";"),
            type_key=paste(unique(type_key), collapse=";"),
            mapped_sex=paste(unique(sort(mapped_sex)), collapse=";"))  %>%
  mutate(mapped_sex=ifelse(mapped_sex=="female;male", "mixed", mapped_sex))

sg_tab_mapped <- single_row %>% bind_rows(mult_row_condensed)
stopifnot(nrow(sg_tab_mapped)==length(sg_tab_mapped$sample_acc))

sg_tab_mapped %>% write_csv("data/01_metadata/mapped_sl_all.csv")


# --- 2. cell line ---- #

# A) using keys
cell_data <- all_attrib_clean %>% 
  filter((str_detect(key_clean, "\\bcell") |
            str_detect(value_clean, "\\bcell")) & 
           !str_detect(key_clean, "tissue")) 

# separately for mouse/human??

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
sample_src_tab <- xenografts %>% select(sample_acc) %>% mutate(xenograft=TRUE) %>%
  full_join(stem_cells %>% select(sample_acc) %>% mutate(stem_cell=TRUE)) %>%
  full_join(primary_cells %>% select(sample_acc) %>% mutate(primary_cell=TRUE)) %>%
  full_join(tissue %>% select(sample_acc) %>% mutate(tissue=TRUE)) %>%
  full_join(cancer %>% select(sample_acc) %>% mutate(cancer=TRUE)) %>%
  full_join(cell_line %>% bind_rows(cell_line2) %>% select(sample_acc) %>% mutate(cell_line=TRUE)) %>%
  full_join(hc_cells %>% select(sample_acc) %>% mutate(cells=TRUE))

sample_src_tab2 <- sample_src_tab %>%
  mutate(across(-sample_acc, ~replace_na(.,FALSE))) %>%
  distinct()
stopifnot(length(unique(sample_src_tab2$sample_acc))==nrow(sample_src_tab2))

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
  mutate(source_type=ifelse(sample_acc %in% blood$sample_acc | 
                               sample_acc %in% tiss_from_cell_line$sample_acc, 
                            "tissue", source_type))

table(src_tab$source_type)

src_tab %>% filter(tissue, cell_line) %>% nrow()

run_to_sample <- read_csv("data/sra_run_to_sample.csv")
compare_src <- src_tab %>% 
  select(sample_acc, source_type) %>%
  full_join(run_to_sample, by=c("sample_acc"="run")) %>%
  mutate(source_type=replace_na(source_type, "other")) %>%
  left_join(sample_type_df %>% rename(metasra_type=sample_type), 
            by=c("samples"="sample_accession")) %>%
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

