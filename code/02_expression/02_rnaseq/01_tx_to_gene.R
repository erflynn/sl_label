# ----  convert transcripts to genes ----- #
# updated version of quant_to_csv that uses tximport to summarize

library('tximport')
library('data.table')
require('tidyverse')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]
idx <- as.numeric(args[2])

extractChunk <- function(my.list, idx, SIZE.CHUNK=50){
  NUM.CHUNKS <- ceiling(length(my.list)/SIZE.CHUNK)
  end_idx <- ifelse((NUM.CHUNKS-1) == idx ,length(my.list), (idx+1)*SIZE.CHUNK)
  return(my.list[(idx*SIZE.CHUNK+1):end_idx])
}

tx2gene <- fread(sprintf("data/rnaseq/transcriptome_indices/%s/long/genes_to_transcripts.txt", prefix),
                              header=FALSE, data.table = FALSE)
colnames(tx2gene) <- c("GENEID", "TXNAME")
tx2gene <- tx2gene %>% select(TXNAME, GENEID)

list.studies <- list.files(sprintf("data/rnaseq/%s/00_infiles/", prefix), pattern="*RP*")
small.list <- extractChunk(list.studies, idx, SIZE.CHUNK=500)


extractGeneCounts <- function(study_id){
  print(study_id)
  files <- list.files(path=sprintf("data/rnaseq/%s/00_infiles/%s/", prefix, study_id), pattern="^[E|S|D]RR[0-9]+")
  names(files) <- str_extract(files, "^[E|S|D]RR[0-9]+")
  if (length(files)==0){
    return(NA)
  }
  files2 <- lapply(files, function(x) sprintf("data/rnaseq/%s/00_infiles/%s/%s", prefix, study_id, x))
  print(head(files2))
  # chunk files if needed 
  num_chunks <- ceiling(length(files2)/50)
  txi_df <- do.call(cbind, lapply(1:num_chunks, function(i){
    list_samples <- extractChunk(files2, i-1)
    txi=tryCatch({
      tximport(list_samples, type = "salmon", tx2gene = tx2gene)
    }, error=function(e){
      file_lengths <- sapply(list_samples, function(x) read.table(pipe(sprintf("wc -l '%s'", x)))[[1]])
      list_samples2 <- list_samples[file_lengths==max(file_lengths)]
      tximport(list_samples2, type = "salmon", tx2gene = tx2gene)
    })

    txi$counts
  }))

  txi_df %>% write_csv(sprintf("data/rnaseq/%s/01_count_mat/%s.csv", prefix, study_id))
}

lapply(small.list, extractGeneCounts)