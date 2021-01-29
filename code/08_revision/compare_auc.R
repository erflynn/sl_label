
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

hm_res <- read_csv("data/hm_alt_class.csv")
hm_res$X1 <- NULL
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
