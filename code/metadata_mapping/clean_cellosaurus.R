require('tidyverse')
require('data.table')


cell_info_df <- fread("../labeling/geo2drug/data/00_db_data/cellosaurus_df.txt", data.table=FALSE)
cell_lab <- cell_info_df %>% separate_rows(synonyms, sep="\\|") %>% select(synonyms, accession, cl)
cell_lab2 <- data.frame(apply(cell_lab, c(1,2) , function(x) str_trim(tolower(x))))
cell_lab_syn <- cell_lab2 %>% select(accession, synonyms) %>% filter(synonyms != "") 
cell_lab_name <- cell_lab2 %>% select(accession, cl) %>% distinct()
cell_df <- rbind(dplyr::rename(cell_lab_syn, cl=synonyms), cell_lab_name) 
cell_df$nwords <- sapply(cell_df$cl, function(x) length(strsplit(x, " ")[[1]]))
cell_df$numchar <- sapply(cell_df$cl, nchar)
cell_df %>% write_csv("data/00_db_data/cell_syn_df.csv")
