# goal is to create a dataset that works well lec3
# topics:
# - window transformations
# lag, lead, slide_vec, cumsum
# - repeating opperations on columns
# - tidying data (wide/long)
# - relational

require('tidyverse')

samples_cl <- read_csv("data/cl_sex_mapped.csv")
study_dates <- read_csv("data/study_dates.csv")

cl <- samples_cl %>% 
  filter(organism == "human" & !is.na(cl_acc)) %>%
  select(sample_acc, cl_acc, cl_name, study_acc) 


cl2 <- cl %>% 
  separate_rows(study_acc, sep=";") %>% 
  left_join(study_dates, by="study_acc")

cl3 <- cl2 %>% 
  group_by(sample_acc) %>% 
  mutate(sample_date=min(date)) %>%
  ungroup() %>%
  select(-study_acc) %>%
  unique() %>%
  filter(!is.na(date))

# --- data frame 1 --- #
cell_line_samples_by_year <- cl3 %>% 
  mutate(year=lubridate::year(date)) %>%
  group_by(year) %>%
  summarize(num_cell_line=n(), 
            num_hela=sum(cl_acc=="cvcl_0030"),
            num_hepg2=sum(cl_acc=="cvcl_0027")) 

# problems
#scores %>%
#  mutate(daily_improvement = score - lag(score))
cl_counts <- cell_line_samples_by_year %>%
  select(year, num_cell_line)
cl_counts %>%
  mutate(increase_by_year=num_cell_line-lag(num_cell_line))


#
#scores = tibble(
#  week = c(1,1,2,2,3,3,4,4),
#  weekday = c("M","F","M","F","M","F","M","F"),
#  score = c(72, 71, 80, 87, 94, 82, 99, 98)
#)

hela_hepg2_by_year <- 
  cell_line_samples_by_year %>% 
  select(year, num_hela, num_hepg2) %>%
  rename(hela=num_hela, hepg2=num_hepg2) %>%
  pivot_longer(cols=c("hela", "hepg2"), 
               names_to="cell_line", 
               values_to="num_samples") %>%
  filter(year > 2010)

# EXERCISE 1: lag/lead
#1. scores %>% mutate(improvment=score-lag(score, n=2))
#2. scores %>% group_by(weekday) %>% mutate(improvement=score-lag(score))
#3. scores %>% group_by(weekday) %>% mutate(improvement=score-lag(score, default=score[1]))
hela_hepg2_by_year %>% mutate(increase_by_year=num_samples-lag(num_samples, n=2))
hela_hepg2_by_year %>% group_by(cell_line) %>% mutate(increase_by_year=num_samples-lag(num_samples))
hela_hepg2_by_year %>% group_by(cell_line) %>% mutate(increase_by_year=num_samples-lag(num_samples, default=num_samples[1]))


cl_counts %>%
  mutate(avg_2_yr = slide_vec(num_cell_line, mean, .before=1))

cl_counts %>%
  mutate(total_samples_to_date=cumsum(num_cell_line))

cl_counts %>%
  mutate(total_samples_to_date = slide_vec(num_cell_line, sum, .before=1))

cl_counts %>%
  mutate(avg_samples_per_yr = slide_vec(num_cell_line, mean, .before=1))


library(weatherData)
sfo = SFO2013 %>%
  transmute( # mutate, but drops all other columns
    time = lubridate::ymd_hms(Time), 
    temp = Temperature
  ) %>%
  group_by(day = lubridate::date(time)) %>%
  summarize(max_temp = max(temp))

# add a column that has the total number of cell lines in the last twelve months
month_year_counts <- cl3 %>%
  mutate(year=lubridate::year(date),
         month=lubridate::month(date)) %>%
  mutate(month_year=sprintf("%s-%s", month, year)) %>%
  group_by(month_year) %>%
  summarize(num_cell_line=n()) 

# answer // TODO
#month_year_counts %>% 

# //TODO could make col-wise part this

# EXERCISE: FILTERING OR REPLACING NAs
#cl2 %>% group_by(cl_acc) %>% count() %>% filter(n>100) %>%ungroup() %>% sample_n(20)

#cl3 %>% #filter(cl_acc=="cvcl_1267") %>%
#  mutate(year=lubridate::year(date)) %>%
#  group_by(year) %>%
#  summarize(num_cell_line=n(),
#            num_hl60=sum(cl_acc=="cvcl_1267"))
 
# // TODO: expression dataset?  `gtex_data()` is a good example for this!
gtex = read_tsv('https://raw.githubusercontent.com/alejandroschuler/r4ds-courses/advance-2020/data/gtex.tissue.zscores.advance2020.txt')
gtex %>% filter()

# MESSY DATA
messy_cl_counts <- cl3 %>% 
  mutate(year=lubridate::year(date)) %>%
  filter(year>=2010) %>%
  group_by(year) %>%
  summarize(hela=sum(str_detect(cl_acc,"cvcl_0030")),
            mcf7=sum(str_detect(cl_acc,"cvcl_0031")),
            hepg2=sum(str_detect(cl_acc, "cvcl_0027")),
            mdamb231 = sum(str_detect(cl_acc, "cvcl_0062")),
            heprg = sum(str_detect(cl_acc, "cvcl_9720"))) %>%
  pivot_longer(cols=hela:heprg, names_to="cell_line", values_to="num_samples") %>%
  pivot_wider(id_cols=c("cell_line"), names_from="year", values_from="num_samples")

tidy_cl_counts <- messy_cl_counts %>%
  pivot_longer(cols=`2010`:`2019`, names_to="year", values_to="num_samples") #%>%
  #mutate(year=str_replace(year, "year", ""),
  #       cell_line=str_replace(cell_line, "num_", "")) %>%
  #select(year, everything()) 
tidy_cl_counts %>%
  ggplot() +
  geom_bar(aes(x=year, y=num_samples, fill=cell_line), stat='identity')

# // TODO: exercise, cleaning GTEx
# // TODO: exercise, creating a table

# Multi-pivoting 
cl_ex <- cl3 %>%
  mutate(year=lubridate::year(date),
         month=lubridate::month(date)) %>%
  mutate(month_year=sprintf("%s-%s", month, year)) %>%
  group_by(month_year) %>%
  summarize(hela=sum(str_detect(cl_acc,"cvcl_0030")),
            mcf7=sum(str_detect(cl_acc,"cvcl_0031")),
            hepg2=sum(str_detect(cl_acc, "cvcl_0027")),
            mdamb231 = sum(str_detect(cl_acc, "cvcl_0062")),
            heprg = sum(str_detect(cl_acc, "cvcl_9720"))) %>%
  ungroup() %>% 
  pivot_longer(cols=hela:heprg, names_to="cell_line", values_to="num_samples") %>%
  filter(str_detect(month_year, "2016|2017"),
                             str_detect(month_year, "9|10|11|12")) %>%
  mutate(month_year=str_replace_all(month_year, "9", "Sept")) %>%
  mutate(month_year=str_replace_all(month_year, "10", "Oct")) %>%
  mutate(month_year=str_replace_all(month_year, "11", "Nov")) %>%
  mutate(month_year=str_replace_all(month_year, "12", "Dec")) %>%
  pivot_wider(id_cols="cell_line", 
              names_from="month_year", 
              values_from="num_samples")


cl_ex %>%
  pivot_longer(cols=contains("-201"), # selects columns that contain this
             names_pattern = "(\\D+)-(\\d+)", # a "regular expression"- we'll learn about these later
             names_to = c(".value", "year")
)


# RELATIONAL DATA EXAMPLE

# -- using GTEx -- #
subject_data <- read_tsv("~/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
sample_data <- read.delim("~/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", stringsAsFactors = FALSE) %>%
  as_tibble()


subj_data2 <- 
  subject_data %>%
  mutate(SEX=case_when(
    SEX==2 ~ "female",
    SEX==1 ~ "male"),
    DEATH=case_when(
         DTHHRDY==0 ~ "ventilator",
         DTHHRDY==4 ~ "terminal illness",
         DTHHRDY==3 ~ "illness",
         DTHHRDY==2 ~ "sudden but natural causes",
         DTHHRDY==1 ~ "accident"
         )) %>%
  select(-DTHHRDY) %>%
  rename(subject_id=SUBJID,
         sex=SEX,
         age=AGE,
         death=DEATH)

sample_data2 <- 
  sample_data %>%
  select(SAMPID, 
         SMNABTCHD,
         #SMGEBTCHD,
         SMCENTER,
         SMTS,
         #SMTSD,
         SMRIN) %>%
  rename(sample_id=SAMPID,
         date_extracted=SMNABTCHD,
         center_id=SMCENTER,
         tissue=SMTS,
         rin_score=SMRIN)
sample_data2$subject_id=sapply(sample_data2$sample_id, function(x) {
  y <- strsplit(x,"-")[[1]];
  sprintf("%s-%s", y[1], y[2])
  }
  )
sample_data2$sample_id=sapply(sample_data2$sample_id, function(x) {
  y <- strsplit(x,"-")[[1]];
  paste(y[3:length(y)], collapse="-")
}
)

sample_data2 <- sample_data2 %>%
  select(subject_id, sample_id, everything()) 
sample_data2 %>%
  mutate(across(everything(),
                ~na_if(., "")))

stopifnot(nrow(sample_data2)==length(unique(sample_data2$sample_id)))


# SLIDE
sample_data2 %>% inner_join(subj_data2, by="subject_id")
# we get a warning, fix this by removing labels!

# why does this fail?
sample_data2 %>% inner_join(subj_data2, by="center_id")


# When keys have different names in different dataframes, the syntax to join is:
gtex_data %>% inner_join(subj_data2, by=c("Ind"="subject_id"))

# EXERCISE:
# find the distinct // TODO

# Join on multiple columns!

# EXERCISE:
# find non-existent planes // TODO



# FIX MISSING DATA --> make it all NAs
#sample_data3 <- sample_data
#  mutate_all(ifelse(.=="",NA, .))
sample_data3 %>% 
  mutate(across(na_if(., ""))
  

