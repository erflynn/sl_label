# process_cell_df.R
# E Flynn
#
# process the cellosaurus data into a data frame
#  filters to only include human, mouse, and rat cell lines!

require('rjson')
require('tidyverse')

cellosaurus <- fromJSON(file="data/00_db_data/cellosaurus_v3.json") 

cell_info_df <- do.call(rbind, 
                        lapply(cellosaurus, function(x) 
                          data.frame(lapply(x, function(y) paste(y, collapse="|")), 
                                     stringsAsFactors=FALSE)))
# note this fails at parsing allele data
cell_info_df$primary_accession <- sapply(cell_info_df$accession, function(x) strsplit(x, "\\|")[[1]][[1]])

cell_info_df$cl <- names(cellosaurus)

# filter for only human, mouse, rat
cell_info_df <- separate_rows(cell_info_df,species, sep="\\|")
cell_info_df2 <- filter(cell_info_df, species %in% c("Homo sapiens", "Mus musculus",  "Rattus norviegicus"))

write.csv(cell_info_df2, "data/00_db_data/cellosaurus_df_v3.txt", row.names=FALSE)


# --- deal with allele data separately --- #
my_dat <- data.frame(do.call(rbind, cellosaurus))
lens <- lapply(my_dat$alleles, length)
cl_names <- rownames(my_dat)[lens!=0]
allele_df <- do.call(rbind, lapply(cl_names, function(cl_name){
  cl_acc = my_dat$accession[[cl_name]]
  l <- my_dat$alleles[[cl_name]]
  df <- do.call(rbind, lapply(names(l), function(allele){
    return(list(
      "cl_acc"=cl_acc[[1]],
      "cl_name" = cl_name,
      "allele" = allele,
      "src" = paste(l[[allele]][["src"]], collapse=";"),
      "num_src"=length(l[[allele]][["src"]])
    ))
  }))
  return(df)
}))
allele_df2 <- data.frame(allele_df, stringsAsFactors = FALSE)
allele_df3 <- apply(allele_df2 , c(1,2), unlist) %>% as_tibble()
allele_freq <- allele_df3 %>% filter(src !="") %>% select(-src) %>% 
  mutate(num_src=as.numeric(num_src)) %>%
  pivot_wider(names_from=allele, values_from=num_src, values_fill=0) %>%
  mutate(num_srcs=`X,Y`+`X`+`Not_detected`+`Y`)
allele_freq %>% write_csv("data/00_db_data/cellosaurus_allele_freq.csv")
# ------ #

