

# code for grabbing RNA-seq data
require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}


studyMat <- function(study_id){
print(study_id)
   out.path <- sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id)
   if (file.exists(out.path) &  file.info(out.path)$size!=0){
       print("file already generated")
       return(NA)
  }

  sample_list <- mapping %>% filter(study_acc==study_id & !is.na(sample_acc) & !is.na(study_acc))
  sample_dat <- lapply(sample_list$sample_acc, function(sample_id){
    sample.path <- sprintf("data/rnaseq/%s/%s/%s_quant.sf", 
                           prefix, study_id, sample_id);
    print(sample_id);
    if (file.exists(sample.path)){
      my_dat <- read_tsv(sample.path) 
      if (!"Name" %in% colnames(my_dat)){
      	 print(sprintf("corrupted file %s", sample.path))
	 return(NA)
      	 #my_dat <- read.delim(sample.path, skip=1)
	 #colnames(my_dat) <- c("Name", "Length","EffectiveLength","TPM", "NumReads")
      }
      quant_sf <- my_dat %>% dplyr::select(Name, TPM);
      colnames(quant_sf) <- c("gene_name", sample_id);
      
      #if (sample_id==sample_list$sample_acc[[1]]){
      #  return(quant_sf)
      #}
      
      return(quant_sf)
    } else {
      return(NA)
    }
  })
  sample_dat2 <- sample_dat[!is.na(sample_dat)]
  if (all(is.na(sample_dat))){
     return(NA)
  }
  row_lengths <- sapply(sample_dat2, nrow)
  if (all(row_lengths==row_lengths[1])){
     study_d1 <- data.frame(do.call(cbind, lapply(sample_dat2, function(df) df %>% select(-gene_name))))
     study_df <- cbind("gene_name"=sample_dat2[[1]]$gene_name, study_d1)
  } else {
    print(sprintf("Row lengths do not match for %s: %s", prefix, study_id))
    study_df <- sample_dat2 %>% reduce(full_join, by="gene_name")
  }
  study_df %>% 
    write_csv(out.path)
}


# read in the mapping file
mapping <- fread(sprintf("data/%s/02_sample_lists/%s_rnaseq_exp_to_sample.csv", prefix, prefix), 
                 data.table=FALSE) %>% filter(!is.na(study_acc))
#list.studies <- unique(mapping$study_acc)
list.studies <- list.files(sprintf("data/rnaseq/%s/", prefix), pattern="*RP*")
print(length(list.studies))
small.list <- extractChunk(list.studies, idx, SIZE.CHUNK=500)
print(small.list)
lapply(small.list, studyMat)

