

# code for grabbing RNA-seq data
require('tidyverse')
require('data.table')
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- args[2]

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}


studyMat <- function(study_id){
  sample_list <- mapping2 %>% filter(study_acc==study_id)
  sample_dat <- lapply(sample_list$sample_acc, function(sample_id){
    sample.path <- sprintf("data/rnaseq/%s/%s/%s_quant.sf", 
                           prefix, study_id, sample_id);
    print(sample_id);
    if (file.exists(sample.path)){
      quant_sf <- read_tsv(sample.path) %>% dplyr::select(Name, TPM);
      colnames(quant_sf) <- c("gene_name", sample_id);
      
      if (sample_id==sample_list$sample_acc[[1]]){
        return(quant_sf)
      }
      
      return(quant_sf[,2])
    } else {
      return(NA)
    }
  })
  
  study_df <- data.frame(do.call(cbind, sample_dat[!is.na(sample_dat)]))
  study_df %>% 
    write_csv(sprintf("data/rnaseq/%s/01_study_mat/%s.csv", prefix, study_id))
}


# read in the mapping file
mapping <- fread(sprintf("data/%s/02_sample_lists/%s_rnaseq_exp_to_sample.csv", prefix, prefix), 
                 data.table=FALSE)
list.studies <- unique(mapping$study_acc)
small.list <- extractChunk(list.studies, idx, 100)
lapply(small.list, studyMat)

