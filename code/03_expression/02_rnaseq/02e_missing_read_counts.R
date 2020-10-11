

require('tidyverse')
require('data.table')

args <- commandArgs(trailingOnly=TRUE)
prefix <- args[1]

my_df <- read_csv("data/missing_read_counts.csv") %>% 
  filter(organism==prefix) %>%
  select(-organism, -num_reads)

getNR <- function(study_id, sample_id){
  sample.path <- sprintf("data/rnaseq/%s/00_infiles/%s/%s_quant.sf", 
                         prefix, study_id, sample_id);
  my_dat <- fread(sample.path, data.table=FALSE) 
  sum(my_dat$NumReads)
}
 
num_reads <- apply(head(my_df), 1, function(x) num_reads=getNR(x[1], x[2]))
my_df$num_reads <- num_reads
write_csv(my_df, sprintf("data/%s_read_counts_m.csv", prefix))

