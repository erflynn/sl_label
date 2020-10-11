# goal:
#  - source type labels
#  - compare to MetaSRA
#  - main goal: use to filter cell line vs non-cell line for creating the test set

require('tidyverse')
load("data/01_metadata/all_sample_attrib_clean.RData") # --> all_attrib_clean

all_attrib_clean2 <- all_attrib_clean %>% 
  dplyr::select(-key, -value) %>% 
  dplyr::rename(key=key_clean, value=value_clean)
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
           str_detect(value, "culture|passage|\\bcell|atcc|huvec"))
cell_line <- all_attrib_clean2 %>% 
  filter(str_detect(key, "cell"),
         str_detect(key, "line"))

cell_line2 <- all_attrib_clean2 %>% 
  anti_join(cell_line) %>%
  filter(str_detect(value, "\\bcell"),
         str_detect(value, "\\bline"))


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
  full_join(cell_line %>% bind_rows(cell_line2) %>% select(sample_acc) %>% 
              mutate(cell_line=TRUE)) %>%
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
src_tab %>% write_csv("data/updated_sample_src.csv")
