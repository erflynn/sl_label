# code for grabbing xy genes 
require('tidyverse')
require('data.table')


args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

load(sprintf("data/rnaseq/%s/03_model_in/expr_sl_in.RData", prefix)) # --> all_df

gene_names <- all_df$gene_name
rm(all_df)

# load the transcript to gene mapping
map_l <- fread(sprintf("data/rnaseq/transcriptome_indices/%s/long/genes_to_transcripts.txt", prefix), 
               header=FALSE, data.table = FALSE)

colnames(map_l) <- c("gene", "transcript")


load(sprintf("../sex_labeling/geo_pipeline/gpl_ref/%s_gene_map.RData", prefix, prefix))
xy_genes <- ref_dat %>% 
  filter(chromosome_name %in% c("X", "Y")) %>% 
  select(chromosome_name, ensembl_gene_id) %>% 
  unique()

xy_transcripts <- inner_join(xy_genes, map_l, by=c("ensembl_gene_id"="gene"))

overlap.transcripts <- xy_transcripts %>% filter(transcript %in% gene_names)
write_csv(overlap.transcripts, sprintf("data/rnaseq/%s/03_model_in/xy_genes_rnaseq.csv", prefix))
