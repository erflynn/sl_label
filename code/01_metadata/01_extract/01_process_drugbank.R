# process_drugbank.R
# E Flynn
# 3/27/2019
# 
# Code for processing drugbank data --> data frame (condensed format) + vocab (long format)
#   - filters out: nutraceuticals, allergens
# TODO: consider filtering more


require('rjson')
require('tidyverse')
source('code/utils/drug_count_utils.R')


# load drugbank data and convert it into a data frame with one row per drug
drugbank <- fromJSON(file="data/00_db_data/drugbank_info.json") 
list.drugs <- lapply(drugbank, function(x) unique(unlist(tolower(c(x$synonyms, x$name)))))
names(list.drugs) <- lapply(drugbank, function(y) y$dbID)

# collapse the categories differently?


drugbank_df <- do.call(rbind,
                       lapply(drugbank, function(x) data.frame(
                         lapply(x, function(y) paste(unlist(y), collapse="|")), stringsAsFactors=FALSE)))
# escape quotes
drugbank_df$name <- sapply(drugbank_df$name, function(x) gsub( '\\"', "\\\'\'", x))
drugbank_df$synonyms <- sapply(drugbank_df$synonyms, function(x) gsub( '\\"', "\\\'\'", x))
write.table(drugbank_df, file="data/00_db_data/drugbank_parsed.txt", sep="\t", row.names=FALSE)

# remove data I don't want! nutraceuticals + allergens
drugbank_groups <- separate_rows(drugbank_df, group, sep="\\|")
drugbank_categories <- separate_rows(drugbank_df, categories, sep="\\|")
cat_table <- table(drugbank_categories$categories)
drugbank_categories$allergen <- drugbank_categories$categories=="Non-Standardized Food Allergenic Extract"
#drugbank_categories$aa <- drugbank_categories$categories %in% c("Amino Acids, Basic", "Amino Acids, Essential")

drugbank_aa_allergen <- drugbank_categories %>% group_by(dbID, name) %>% summarize(discard=any(allergen))
discard_dbs <- drugbank_aa_allergen$dbID[drugbank_aa_allergen$discard]
drugbank_no_nutra <- filter(drugbank_groups, group!="nutraceutical" & !(dbID %in% discard_dbs))

# --- CONSTRUCT DRUGBANK VOCABULARY --- #
drug_name_syn <- createNameSynVocab(drugbank_no_nutra, "dbID")

write.table(drug_name_syn, file="data/00_db_data/drugbank_vocab_no_nutra.txt", sep="\t", row.names=FALSE)
