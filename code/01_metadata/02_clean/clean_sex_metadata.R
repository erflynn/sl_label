# clean_sex_metadata.R
# E Flynn
#
# Code to grab sex annotations that don't fall into "male" or "female"
# These are then manually mapped to male, female, mixed, or unknown

human_metadata <- read.csv("data/01_metadata/human_metadata.csv")
mouse_metadata <- read.csv("data/01_metadata/mouse_metadata.csv")
rat_metadata <- read.csv("data/01_metadata/rat_metadata.csv")

# human_metadata %>% filter(!is.na(sex)) %>% select(sex) %>% unique() %>% write_csv("human_alt_sex_annot.csv")
# mouse_metadata %>% filter(!is.na(sex)) %>% select(sex) %>% unique() %>% write_csv("mouse_alt_sex_annot.csv")
mouse_alt <- read_csv("mouse_alt_sex_annot.csv")
human_alt <- read_csv("human_alt_sex_annot.csv")

