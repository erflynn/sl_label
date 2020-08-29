# ran to check differences, does not need to repeat

require("tidyverse")

# --- get arrayexpress list --- #
array_exp <- comb_metadata %>% filter(!str_detect(sample_acc, "GSM|SRR|ERR|DRR"))
array_exp_studies <- array_exp %>% distinct(study_acc)
array_exp_studies %>% write_tsv("data/array_express_studies.tsv", col_names=FALSE)


list.studies <- list.files(sprintf("data/rnaseq/%s/", prefix), pattern="*RP*")

print(length(list.studies))
list.studies <- list.studies[list.studies != "README.md"]
small.list <- extractChunk(list.studies, idx, SIZE.CHUNK=500)

meta_res <- lapply(list.studies, function(study_id){
  my.f <- sprintf("data/rnaseq/%s/%s/metadata_%s.tsv", prefix, study_id, study_id)
  if (!file.exists(my.f)){
    return(NA)
  }
  meta_ind <- read_tsv(my.f)
  return(meta_ind)
})
meta_df <- do.call(rbind, meta_res[!is.na(meta_res)])

# load all the metadata and compare to the main metadata file

metadata <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix))
mapping <- read.csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample.csv", prefix))

setdiff(metadata$acc, meta_df$refinebio_accession_code)
# ok so these are the same -- but the mapping file does overdo things, but that's ok.

