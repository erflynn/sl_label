# Compare sample source type data 

require('tidyverse')

src_tab <- read_csv("data/updated_sample_src.csv")
run_to_sample <- read_csv("data/sra_run_to_sample.csv")
load("data/10_comparison/metasra_data.RData") # --> sample_type_df
# //TODO: check

compare_src <- src_tab %>% 
  dplyr::select(sample_acc, source_type) %>%
  full_join(run_to_sample, by=c("sample_acc"="run")) %>%
  mutate(source_type=replace_na(source_type, "other")) %>%
  left_join(sample_type_df %>% dplyr::rename(metasra_type=sample_type), 
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
  dplyr::count() %>%
  pivot_wider(names_from="source_type", values_from="n", values_fill=0)
compare_src %>% 
  filter(!metasra_type %in% c( "other" , "in vitro differentiated cells") &
           !source_type %in% c("other", "xenograft", "cell")) %>%
  mutate(match=(metasra_type==source_type)) %>%
  group_by(match) %>% dplyr::count() 

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
  group_by(metasra_type, source_type) %>% dplyr::count() %>%
  pivot_wider(names_from=source_type, values_from=n) 
simple_comp %>% 
  filter(metasra_type %in% c("cell", "tissue", "other"),
         source_type %in% c("cell", "tissue", "other")) %>%
  mutate(across(c(metasra_type, source_type), 
                ~factor(.,levels=c("cell", "tissue", "other")))) %>%
  mutate(tot=n()) %>%
  group_by(metasra_type, source_type, tot) %>% dplyr::count() %>%
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

library('rjson')
metasra_sample_type =fromJSON(file="data/00_db_data/metasra_training_set_sample_type.json")

extractAccSample <- function(injson){
  accs <- sapply(injson, function(x) x$sample_accession)
  sample_types <- sapply(injson, function(x) x$sample_type)
  df_dat <- tibble("sample"=accs, "source_type"=sample_types) %>% 
    left_join(run_to_sample, by=c("sample"="samples"))
  print(table(is.na(df_dat$run)))
  df_dat2 <- df_dat %>% filter(!is.na(run))
  compare_df <- df_dat2 %>% 
    mutate(source_type=str_replace_all(source_type, "_", " ")) %>%
    left_join(src_tab  %>% select(sample_acc, source_type),
              by=c("run"="sample_acc")) 
  
  return(compare_df)
}
metasra_training <- extractAccSample(metasra_sample_type)

# more metasra training
metasra_sample_type_test1 =fromJSON(file="data/00_db_data/metasra_test_nonenriched_sample_type.json")
#metasra_sample_type_test2 =fromJSON(file="data/00_db_data/metasra_test_enriched_sample_type.json")
metasra_test1 <- extractAccSample(metasra_sample_type_test1)
#metasra_test2 <- extractAccSample(metasra_sample_type_test2)
table(metasra_test1$source_type.x, metasra_test1$source_type.y)
#table(metasra_test2$source_type.x, metasra_test2$source_type.y)
#>> not including enriched b/c it is identical to non-enriched


table(metasra_training$source_type.x==metasra_training$source_type.y) # 400

cl_or_tiss <- metasra_training %>% filter(source_type.y %in% c("cell line", "cell", "tissue", "primary cells"))

cl_or_tiss %>% filter( source_type.y=="primary cells")
all_attrib_clean2 %>% filter(sample_acc=="SRR309267")
# this seems HAZY af

cl_or_tiss <- metasra_training %>% 
  filter(source_type.y %in% c("cell line", "tissue"))

table(cl_or_tiss$source_type.x, cl_or_tiss$source_type.y)
all_attrib_clean2 %>% 
  semi_join(cl_or_tiss %>% 
              filter(source_type.x=="cell line" & source_type.y=="tissue"), 
            by=c("sample_acc"="run")) %>% View()

