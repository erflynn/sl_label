# Set up code for running HC pipeline
#

require('tidyverse')
by_study <- read_csv("data/study_sex_lab.csv")
comp_metadat


study_input_hc <- by_study %>% 
  filter(label_type=="metadata",
         study_sex=="mixed sex" ,
         num_present >=10 & num_f >= 4 & num_m >= 4)  

study_input_hc %>%
  group_by(organism, data_type) %>%
  count()
