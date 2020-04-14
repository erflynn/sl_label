# process_cell_df.R
# E Flynn
#
# process the cellosaurus data into a data frame
#  filters to only include human, mouse, and rat cell lines!

require('rjson')
require('tidyverse')

cellosaurus <- fromJSON(file="data/00_db_data/cellosaurus_v2.json") 

cell_info_df <- do.call(rbind, 
                        lapply(cellosaurus, function(x) 
                          data.frame(lapply(x, function(y) paste(y, collapse="|")), 
                                     stringsAsFactors=FALSE)))
# CHANGE THIS

cell_info_df$primary_accession <- sapply(cell_info_df$accession, function(x) strsplit(x, "\\|")[[1]][[1]])

cell_info_df$cl <- names(cellosaurus)

# filter for only human, mouse, rat
cell_info_df <- separate_rows(cell_info_df,species, sep="\\|")
cell_info_df2 <- filter(cell_info_df, species %in% c("Homo sapiens", "Mus musculus",  "Rattus norviegicus"))

write.csv(cell_info_df2, "data/00_db_data/cellosaurus_df_v2.txt", row.names=FALSE)
