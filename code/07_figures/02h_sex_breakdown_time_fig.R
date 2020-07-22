# 02h_sex_breakdown_time_fig.R
# E Flynn
# 07/22/2020
#
# Make a figure showing sex breakdown over time.
# Relies on extracting array express (02b-d) and RNA-seq (02e-g) data first.
#  We do not need to extract GEO data because it is already present
# 
# First construct these tables:
# 1. sample | organism | data_type | sex | studies | sample_date
# 2. study | study_date
#
# Note: skipping sample-level goal, and instead going to focus on just the study/study-date
# TODOs (remaining):
# - where do mostly male + mostly female go? for now ignoring
# - add in more tests to make sure correct
# - make legend prettier, consider putting together?
# - calculate a stat about change over time
# - aggregate panel for RNA-seq + microarray?

require('tidyverse')
require('lubridate')
require('zoo')


comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv")

h_metadata <- read.csv("data/01_metadata/human_experiment_metadata.csv", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  filter(study_acc %in% list_studies) %>%
  select(study_acc, date) 

m_metadata <- read.csv("data/01_metadata/mouse_experiment_metadata.csv", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  filter(study_acc %in% list_studies) %>%
  select(study_acc, date)

rb_micro_meta <- h_metadata %>% bind_rows(m_metadata) %>% mutate(date=ymd_hms(date))



ae_dates <- read_csv("data/dates/ae_date_info.csv") %>% 
  select(-diff_yr, -update_date) %>%
  rename(study_acc=accession, date=release_date) %>%
  unique()

sra_dates <- read_csv("data/dates/sra_study_ena.csv") %>%
  rename(study_acc=accession) %>%
  unique()

study_dates2 <- rb_micro_meta %>% 
  mutate(date = ymd(date)) %>% 
  unique() %>%
  filter(!is.na(date)) %>% 
  bind_rows(ae_dates) %>% 
  unique() %>%
  bind_rows(sra_dates) %>%
  unique()

stopifnot(length(unique(study_dates2$study_acc))==nrow(study_dates2))
# check AE filled in missing
stopifnot((length(unique(study_dates2$study_acc))-nrow(sra_dates))==length(unique(rb_micro_meta$study_acc)))

study_sl <- read_csv("data/study_sex_lab.csv")

sl_study2 <- study_sl %>% left_join(study_dates2, by="study_acc") %>%
  mutate(study_sex=factor(study_sex, levels=c("female only", "mostly female", 
                                              "mixed sex", "mostly male", 
                                              "male only", "unknown"))) %>%
  mutate(data_type=case_when(data_type=="RNA-seq" ~ "rnaseq",
                             TRUE ~ data_type)) %>%
  mutate(label_type=factor(label_type, levels=c("metadata", "expression")))

sl_study2 %>% filter(is.na(date)) %>% select(study_acc) %>% unique() %>% nrow()
# only 18 documented are missing these data


# -- original plot -- #
ggplot(sl_study2, aes(x=date))+
  geom_histogram(aes( fill=study_sex), bins=15)+
  ylab("number of studies")+
  xlab("date submitted")+
  facet_grid(rows=vars(label_type), cols=vars(organism),  scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())  


# -- tile over time and calculate fractions -- #
dat_plus2 <- sl_study2 %>%
  arrange(organism, data_type, label_type, study_sex,date) %>%
  mutate(value=1) %>%
  group_by(organism, data_type, label_type, study_sex) %>%
  filter(!is.na(date)) %>%
  complete(date = full_seq(date, period = 1), fill = list(value = 0)) %>%
  group_by(organism, data_type, label_type, study_sex, month3=floor_date(date, "3 months")) %>%
  summarize(value=sum(value)) %>%
  mutate(tot= rollapplyr(value, width = 4, FUN = sum, partial = TRUE)) %>%
  select(organism,data_type, label_type, study_sex, month3, tot) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(id_cols=c(organism, data_type, label_type, month3),
              names_from=study_sex, values_from=tot, values_fn=list(tot=max)) %>%
  mutate_if(is.numeric, ~ifelse(is.na(.), 0, .)) %>%
  group_by(organism, data_type, label_type, month3) %>%
  mutate(num_studies=(`female only`+`male only` + `mostly male` + 
                        `mostly female` + `mixed sex` + `unknown`)) %>%
  mutate(frac_missing=(`unknown`)/num_studies,
    rat_single_sex=(`female only`+`male only`)/(`mixed sex`+1),
    frac_mixed_sex=(`mixed sex`)/num_studies,
    frac_single_sex=(`female only`+`male only`)/num_studies,
    frac_m_only=(`male only`)/num_studies,
    frac_f_only=(`female only`)/num_studies,
    rat_m_only=`male only`/(`female only`+`male only`+1))

# TODO: add in some sanity checks here!

# grab the frac and total data to make a figure
dat_plus3 <- dat_plus2 %>% 
  select(organism, data_type, label_type, month3, num_studies, frac_missing, 
                     frac_mixed_sex, frac_single_sex, frac_m_only, frac_f_only)
dat_plus_long <- dat_plus3 %>% 
  pivot_longer(cols=num_studies:frac_f_only, names_to="col_desc") %>%
  mutate(col_type=ifelse(str_detect(col_desc, "frac"), 
                         "fraction of studies", 
                         "number of studies")) %>%
  ungroup() %>%
  filter((label_type=="expression" & col_desc!="frac_missing") |
           (label_type=="metadata" & col_desc=="frac_missing")) %>%
  filter(col_desc != "frac_single_sex") %>%
  mutate(col_desc=factor(col_desc, levels=c("frac_f_only", 
                                            "frac_missing",
                                            "frac_mixed_sex", 
                                            "frac_m_only",
                                            "num_studies")))

# -- make figure!! -- #
ggplot(dat_plus_long, aes(x=month3, y=value))+
  geom_line(aes(lty=data_type, col=col_desc))+
  facet_grid(vars(col_type), vars(organism), scales="free")+
  theme_bw()+
  
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.title = element_blank()) +
  xlab("")+
  ylab("")
# // TODO reorder legend, colors? make this two plots?
# // TODO: should there be an aggregate panel?
ggsave("figures/paper_figs/fig2.png" )




# -- extra exploratory figs -- #
ggplot(dat_plus2, aes(x=month3, y=rat_single_sex))+
  geom_line(aes(col=organism))+
  facet_grid(rows=vars(data_type), cols=vars(label_type),  scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  xlab("")+
  ylab("ratio of single-sex to mixed-sex")

ggplot(dat_plus2 %>% filter(label_type=="metadata"), aes(x=month3, y=frac_missing))+
  geom_line(aes(col=data_type))+
  facet_grid(.~(organism),  scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  xlab("")+
  ylab("fraction missing")

ggplot(dat_plus2, aes(x=month3, y=frac_missing))+
  geom_line(aes(col=data_type))+
  facet_grid(rows=vars(organism), cols=vars(label_type),  scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  xlab("")+
  ylab("fraction missing")


ggplot(dat_plus2, aes(x=month3, y=rat_m_only))+
  geom_line(aes(col=organism))+
  facet_grid(rows=vars(data_type), cols=vars(label_type),  scales="free")+
  theme_bw()+
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  xlab("")+
  ylab("ratio of male-only to female-only")


