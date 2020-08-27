
require('tidyverse')

# TODO - how did I get these data?
human_compare <- read_csv("data/human_sl_compare.csv") # 17119
mouse_compare <- read_csv("data/mouse_sl_compare.csv") # 8836

compare_dat <- human_compare %>% bind_rows(mouse_compare)

comb_metadata <- read_csv("data/01_metadata/combined_human_mouse_meta.csv", col_types="cccccccdcd")
microarray <- comb_metadata %>% filter(data_type=="microarray") %>% select(-present, -num_reads)
source_type <- read_csv("data/sample_source_type.csv")


compare_df <- compare_dat %>% 
  rename(sample_acc=gsm) %>% 
  inner_join(microarray %>% select(sample_acc, organism, expr_sex, p_male)) %>% # .. hmmm, we get 15426/25955
  inner_join(source_type %>% select(acc, source_type, cl_line), by=c("sample_acc"="acc")) %>% # not lossy
  select(sample_acc, organism, source_type, cl_line, everything())

# TODO - how are these missing?!
missing_gsms <- setdiff(compare_dat$gsm, microarray$sample_acc)

compare_df2 <- compare_df %>% 
  #filter(!source_type %in% c("named_cl", "unnamed_cl")) %>% 
  filter(!is.na(text_sex)) %>%
  filter(!(is.na(toker_sex) & is.na(massir_sex))) %>%
  mutate(compare_col=case_when(
    is.na(massir_sex) & text_sex==toker_sex ~ 1,
    is.na(massir_sex) & text_sex!=toker_sex ~ -1,
    is.na(toker_sex) & massir_sex==text_sex ~ 1,
    is.na(toker_sex) & massir_sex!=text_sex ~ -1,
    toker_sex==text_sex & massir_sex==text_sex ~ 2,
    toker_sex!=text_sex & massir_sex!=text_sex ~ -2,
    toker_sex == text_sex & massir_sex != text_sex ~ 0,
    toker_sex != text_sex & massir_sex == text_sex ~ 0
  ))
  
compare_df2 %>%
  mutate(expr_sex=ifelse(p_male < 0.7 & p_male > 0.3, "unknown", expr_sex)) 

# how many of the matching ones do we call correctly? at each threshold?
matching <- compare_df2 %>% filter(compare_col %in% c(1,2))

summarizeAcc <- function(df, x){
  df %>%
    mutate(threshold=x) %>%
    select(sample_acc, organism, source_type, text_sex, expr_sex, p_male, compare_col, threshold) %>%
    mutate(expr_sex=ifelse(p_male < threshold & p_male > 1-threshold, "unknown", expr_sex)) %>%
    mutate(match=case_when(
      expr_sex=="unknown" ~ 0,
      expr_sex==text_sex ~ 1,
      expr_sex!=text_sex ~ -1)) %>%
    group_by(organism, source_type, threshold,match) %>% 
    count() %>%
    pivot_wider(names_from=match, names_prefix="match_type", values_from=n, values_fill=0) %>%
    mutate(nsamples=sum(`match_type-1`+match_type0+match_type1)) %>%
    mutate(across(`match_type-1`:match_type1, ~./nsamples))
}
  
df <- do.call(rbind, lapply(seq(0.6,0.9, 0.1), function(x) summarizeAcc(matching, x))) %>%
  mutate(frac_unlab=match_type0,
         frac_correct=abs(match_type1-`match_type-1`)/(match_type1+`match_type-1`))
df %>% filter(threshold==0.7) # 90.9% (5.15% unlab) for mouse, 99.5% (3.17% unlab) for human at threshold 0.7

# how many mismatch across all? at each threshold?
mismatch <- compare_df2 %>% filter(compare_col %in% c(-1,-2))

df2 <- do.call(rbind, lapply(seq(0.6,0.9, 0.1), function(x) summarizeAcc(mismatch, x))) 
df3 <- do.call(rbind, lapply(seq(0.6,0.9, 0.1), function(x) summarizeAcc(compare_df2, x))) %>%
  mutate(num_unlab=match_type0*nsamples) %>%
  group_by(organism, threshold) %>%
  summarize(num_unlab=sum(num_unlab), n=sum(nsamples)) %>%
  mutate(frac_unlab=num_unlab/n)
df3
# //TODO - add unlabeled counts
  
mismatch_df2 <- df2 %>% 
  mutate(num_mismatch=`match_type-1`*nsamples) %>% 
  #mutate(frac_unlab=match_type0) %>%
  select(-contains("type"), -contains("frac")) %>%
  rename(n_mismatch_other=nsamples) %>%
  left_join(compare_df %>% group_by(organism, source_type) %>% count()) %>%
  mutate(frac_mismatch=num_mismatch/n)

mismatch_df2 %>% 
  ungroup() %>%
  group_by(organism, threshold) %>%
  summarize(num_mismatch=sum(num_mismatch), n=sum(n)) %>%
  mutate(frac_mismatch=num_mismatch/n)

# summarize = all

mismatch_df2 %>% 
  ungroup() %>%
  mutate(cl_line=(source_type %in% c("unnamed_cl", "named_cl"))) %>%
  group_by(organism, threshold, cl_line) %>%
  summarize(num_mismatch=sum(num_mismatch), n=sum(n)) %>%
  mutate(frac_mismatch=num_mismatch/n) %>%
  arrange(cl_line, organism, threshold)

# summarize = non-cl

