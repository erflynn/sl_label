# count the number of samples and studies, 
# and the overlap for a supplement table


require('tidyverse')
options(stringsAsFactors = FALSE)

prefix <- "human"

compendia <- read.csv(sprintf("data/01_metadata/%s_metadata.csv", prefix)) # 430119
rnaseq <- read.csv(sprintf("data/01_metadata/%s_rnaseq_sample_metadata.csv", prefix)) # 229789
exp_to_sample <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample.csv", prefix))
exp_to_sample_rnaseq <- read_csv(sprintf("data/01_metadata/%s_exp_to_sample_counts.csv", prefix))
exp_to_sample_rnaseq2 <- exp_to_sample_rnaseq %>% filter()

overlap <- intersect(compendia$acc, rnaseq$acc) # 99611 overlapping IDs
microarray <- setdiff(compendia$acc, rnaseq$acc) # 330508
rnaseq_only <- setdiff(rnaseq$acc, compendia$acc) # 130178

num_compendia <- nrow(compendia)
num_overlap <- nrow(overlap)
num_microarray <-


# number of studies

# number of samples

print(nrow(overlap))
print(nrow(microarray))
print(nrow(rnaseq_only))