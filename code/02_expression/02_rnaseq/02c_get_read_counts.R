require('tidyverse')
require('data.table')
args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}


exp_samp <- read_csv(sprintf("data/01_metadata/%s_rnaseq_exp_to_sample2.csv", prefix))
exp_present <- exp_samp %>% filter(present)

my_chunk <- extractChunk(1:nrow(exp_present), idx, 5000)

getNR <- function(study_id, sample_id){
  sample.path <- sprintf("data/rnaseq/%s/00_infiles/%s/%s_quant.sf", 
                         prefix, study_id, sample_id);
  my_dat <- fread(sample.path, data.table=FALSE) 
  sum(my_dat$NumReads)
}
nr_df <- data.frame()
for (i in 0:4){
    my_df <- exp_present[my_chunk[(i*1000+1):min(((i+1)*1000), length(my_chunk))],]
    num_reads <- apply(my_df, 1, function(x) num_reads=getNR(x[1], x[2]))
    my_df$num_reads <- num_reads
    if (i==0){
       nr_df <- my_df 
    } else {
      nr_df <- rbind(nr_df, my_df)
  }
}  

nr_df %>% write_csv(sprintf("data/rnaseq/%s/02_read_counts/%s_read_counts_%s.csv", prefix, prefix, idx))

