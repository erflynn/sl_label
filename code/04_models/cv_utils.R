# helpful utilities for separating studies!


require('tidyverse')

sep_studies <- function(num_studies, nfolds){
  my_l <- c(1:nfolds, nfolds:1)
  if (runif(1, 0, 1) >= 0.5){
    my_l <- c(nfolds:1, 1:nfolds)
  }
  if (num_studies < (2*nfolds)){
    
    return(my_l[1:num_studies])
  } 
  else {
    num_reps <- num_studies %/% (2*nfolds)
    num_rem <- num_studies %% (2*nfolds)
    if (num_rem==0){
      return(rep(my_l, num_reps))
    } else{
      return(c(rep(my_l, num_reps), my_l[1:num_rem]))
    }
  }
}

partition_group_data <- function(df, grp_col="grp", class_col="class", nfolds=2){
  colnames(df)[colnames(df)==grp_col] <- "grp"
  colnames(df)[colnames(df)==class_col] <- "class"
  # get the counts by class in each grp
  study_counts_by_class <- df %>% 
    mutate(grp=as.factor(grp)) %>%
    group_by(grp, class) %>% 
    count() %>% 
    ungroup() %>%
    pivot_wider(names_from=class, values_from=n, names_prefix="num", values_fill=c(n=0)) 
  
  # shuffle and partition
  partitioned_data <- study_counts_by_class %>% 
    group_by_if(is.numeric) %>% 
    sample_n(n()) %>%
    mutate(partition=unlist(sep_studies(n(), nfolds)))  %>%
    ungroup()
  
  # add the sample names back in
  samples_to_grps <- partitioned_data %>% 
    select(grp, partition) %>%
    left_join(df, by=c("grp")) 
  
  return(samples_to_grps)
}


#set.seed(104)
#grp_partitioned_dat <- partition_group_data(df, grp_col ="participant", class_col="diagnosis")
#grp_partitioned_dat %>% group_by(partition, class) %>% count()
