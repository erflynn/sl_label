
library(tidyverse)

auc_plot <- function(df){
  df_new <- df %>% distinct(method, acc, auc) %>%
    mutate(classifier=sprintf("%s (%s", method, round(acc, 3)*100)) %>%
    mutate(`classifier (accuracy)`=paste(classifier, "%)", sep="")) %>%
    select(method, `classifier (accuracy)`)
  ggplot(df %>% left_join(df_new), 
         aes(x=fpr, y=tpr, col=`classifier (accuracy)`))+
    geom_line(alpha=0.6)+
    theme_bw()+
    geom_abline(alpha=0.6, slope=1, lty=2)+
    ylab("TPR")+xlab("FPR")+
    xlim(0,1)+
    ylim(0,1)
}

load_acc <- function(prefix){
  hm_res <- do.call(rbind, lapply(1:6, function(x) 
    read_csv(sprintf("data/%s_alt_class_%s.csv", prefix, x))))
  hm_res$X1 <- NULL
  hm_acc <- hm_res %>% distinct(acc, method) %>%
    arrange(desc(acc))
  return(hm_res)
}



acc_df <- do.call(rbind, 
        list(load_acc("mr") %>% mutate(organism="mouse", data_type="RNA-seq"), 
          load_acc("hr")%>% mutate(organism="human", data_type="RNA-seq"), 
          load_acc("mm")%>% mutate(organism="mouse", data_type="microarray"), 
          load_acc("hm") %>% mutate(organism="human", data_type="microarray")))

acc_df2 <- acc_df %>% distinct(acc, method) %>%
  arrange(desc(acc))
acc_df$method <- factor(acc_df$method, levels=unique(acc_df2$method))
acc_df %>%  
  ggplot(aes(x=method, y=acc))+geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Accuracy")+
  xlab("")+
  facet_grid(data_type ~ organism)
ggsave("figures/revision/ml_methods_comparison.png")

auc_plot(hm_res)
ggsave("figures/revision/hm_alt_models.png")


hr_res <- read_csv("data/hr_alt_class.csv")
hr_res$X1 <- NULL
auc_plot(hr_res)
ggsave("figures/revision/hr_alt_models.png")

mm_res <- read_csv("data/mm_alt_class.csv")
mm_res$X1 <- NULL
auc_plot(mm_res)
ggsave("figures/revision/mm_alt_models.png")

mr_res <- read_csv("data/mr_alt_class.csv")
mr_res$X1 <- NULL
auc_plot(mr_res)
ggsave("figures/revision/mr_alt_models.png")

methods_comp <- mm_res %>% 
  distinct(method, auc, acc) %>%
  mutate(organism="mouse", data_type="microarray") %>%
  bind_rows(mr_res %>% 
              distinct(method, auc, acc) %>%
              mutate(organism="mouse", data_type="rnaseq")) %>%
  bind_rows(hr_res %>% 
              distinct(method, auc, acc) %>%
              mutate(organism="human", data_type="rnaseq")) %>%
  bind_rows(hm_res %>% 
              distinct(method, auc, acc) %>%
              mutate(organism="human", data_type="microarray")) 

methods_comp %>%  arrange(desc(acc)) %>%
  arrange(organism, data_type, desc(acc)) %>%
  group_by(organism, data_type) %>% 
  slice_max(n=3, order_by=acc)
