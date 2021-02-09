# code for extracting the file
# head -1 MUS_MUSCULUS.tsv > list_cols.txt

# get the list of files
# add f_idx,etc

library(vroom)
library(tidyverse)
NUM_PER_FILE <- 20000
lc <- vroom("list_cols.txt", col_names=F)

# rotate
lc_t <- data.frame(t(lc))
colnames(lc_t) <- "sample_acc"
lc_t <- lc_t %>% 
  filter(!is.na(sample_acc))
nrow(lc_t) # 430119 228708

# add f_idx, idx
lc_t <- lc_t %>% 
  mutate(idx=1:n())
lc_t2 <- lc_t %>% mutate(f_idx = idx %/%  NUM_PER_FILE)
lc_t2 %>% write_csv("sample_to_file.csv")


# -- samples to check - #
sf <- read_csv("data/03_expression/microarray/mouse/sample_to_file.csv")
sf2 <- sf %>% mutate(rem_idx=idx %%  NUM_PER_FILE) 
max_f_idx <- max(sf2$f_idx)
last_idx <- tail(sf2$rem_idx, 1)

sf_to_check <- sf2 %>% filter(f_idx %in% c(0, 1, max_f_idx),
               rem_idx %in% c(0, 1, 2, last_idx, NUM_PER_FILE-1)) %>%
  filter(!(rem_idx == last_idx & f_idx %in% c(0,1)))
sf_to_check %>% write_csv("m_samp_to_check.csv")

# check file 0:
#  rem_idx = 1, 2
#  rem_idx = 19999

# file 1:
#  0,1,2, 19999

# last file:
#  0,1,2, last_idx