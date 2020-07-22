# Plot the cell line breakdown and sample switching

require('tidyverse')

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

compendia_cl <- read_csv("data/02_labeled_data/human_compendia_sample_cl.csv")
rnaseq_cl <- read_csv("data/02_labeled_data/human_rnaseq_sample_cl.csv")
head(compendia_cl)
head(rnaseq_cl)
