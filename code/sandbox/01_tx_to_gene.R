# ----  convert transcripts to genes ----- #
library('tximport')
library('data.table')
library('tidyverse')

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
  if (length(files)==0){
    return(NA)
  }
  files2 <- sapply(files, function(x) sprintf("data/rnaseq/%s/00_infiles/%s/%s", prefix, study_id, x))
  names(files2) <- str_extract(files, "^[E|S|D]RR[0-9]+")
  
  # run in chunks of files if there are >50 samples
  num_chunks <- ceiling(length(files2)/50)
  list_txs <- lapply(1:num_chunks, function(i){
    list_samples <- extractChunk(files2, i-1, SIZE.CHUNK=50)
    print(list_samples)
    txi=tryCatch({
      res <- tximport(list_samples, type = "salmon", tx2gene = tx2gene)
      my_df <- data.frame(res$counts)
      my_df$GENEID <- rownames(my_df)
      return(my_df)
    }, error=function(e){ 
      print("error")
      
      # errors often occur when there are multiple different lengths 
      # to fix this, we group by file lengths and load each length separately
      # then we join together (this puts NAs where it is missing)
      
      file_lengths <- sapply(list_samples, function(x) read.table(pipe(sprintf("wc -l '%s'", x)))[[1]])
      
      my_l <- tibble("sample"=names(file_lengths), 
                     "file_length"=file_lengths) %>%
        filter(file_length > 10000)
      txi2 <- lapply(my_l %>% group_split(file_length), function(x)
        {list_samples2 <-list_samples[x$sample   ]
        my_df <-data.frame(tximport(list_samples2, type = "salmon", tx2gene = tx2gene)$counts)
        my_df$GENEID <- rownames(my_df)
        return(my_df)
      })
      df <- txi2 %>% reduce(full_join, by="GENEID") # join by geneid
      return(df)
    })
  })
  
  # put together the chunks with a join
  if (length(list_txs) > 1){
    txi_df <- list_txs %>% reduce(full_join, by="GENEID")
  } else {
    txi_df <- list_txs[[1]]
  }
  
  # write it out
  txi_df %>% select(GENEID, everything()) %>% write_csv(sprintf("data/rnaseq/%s/01_count_mat/%s.csv", prefix, study_id))
}

lapply(small.list, extractGeneCounts)