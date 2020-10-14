# 02d_array_express_df.R
# E Flynn
# 07/20/2020
#
# Grab array express date data
# see 02b, 02c where I downloaded non_gse_studies and parsed XML
# ArrayExpress metadata was grabbed using the following
# https://www.ebi.ac.uk/arrayexpress/help/programmatic_access.html#Experiment_details
# > https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/E-xxxx-nnnnn
# > https://www.ebi.ac.uk/arrayexpress/xml/v3/experiments/E-xxxx-nnnnn/samples
# this is XML, so we had to download and parse


require('tidyverse')
require('GEOmetadb')
require('rjson')
require('lubridate')


non_gse_studies %>% read_csv("data/dates/non_gse.csv")

ae_dates <- fromJSON(file="data/dates/ae_date.json") 
ae_date_df <- do.call(rbind, lapply(ae_dates, function(x) data.frame(x, stringsAsFactors = FALSE)))
ae_date_df$accession <- rownames(ae_date_df)
stopifnot(nrow(ae_date_df)==nrow(non_gse_studies))

ae_date_df2 <- ae_date_df %>% 
  as_tibble() %>%
  select(accession, everything()) %>%
  mutate(release_date=ymd(release_date),
         update_date=ymd(update_date)) %>%
  mutate(diff_date=update_date-release_date) %>%
  mutate(diff_date=as.numeric(diff_date)) %>%
  mutate(diff_yr=diff_date/365)

summary(ae_date_df2)
# no missing data!, mean/median ~5-6
any(is.na(ae_date_df2$release_date))
any(is.na(ae_date_df2$update_date))

ggplot(ae_date_df2, aes(x=diff_yr))+geom_histogram() 
# 22 that are negative... 

ae_date_df2 %>% filter(diff_yr < 0) # :/ 

# WEBSITE uses release_date, I think we should stick w...
ae_date_df2 %>% select(-diff_date) %>% write_csv("data/dates/ae_date_info.csv")
