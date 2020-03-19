
require('tidyverse')
require('ggalluvial')

load("data/sample_level_sex.RData" )
load("data/summary_files.RData")


# ----- look at the distribution of mixed sex sample scores ---- #
ggplot(human2 %>% filter(text_sex!="unknown"), aes(x=pred, group=text_sex))+geom_density(aes(col=text_sex))
ggplot(mouse2 %>% filter(text_sex!="unknown"), aes(x=pred, group=text_sex))+geom_density(aes(col=text_sex))
ggplot(rat2 %>% filter(text_sex!="unknown"), aes(x=pred, group=text_sex))+geom_density(aes(col=text_sex))

# mouse one appears interesting, really unclear what else is going on?

# ----- alluvial plots and fractions ---- #

# sample level
alluvPlotSample <- function(ds){
  flow_freq_counts <- ds %>% 
    ungroup() %>% 
    rename(metadata=text_sex, expression=expr_sex) %>%
    group_by(metadata, expression) %>% 
    mutate(Freq=n()) %>% 
    select(-gsm) %>% 
    unique() %>%
    ungroup() %>%
    mutate(row_id=1:n()) %>%
    gather(key="labeling_method", value="sex", -Freq, -row_id) %>%
    mutate(sex=ifelse(sex=="unknown", "unlabeled", sex)) %>%
    mutate(row_id=as.factor(row_id), 
           labeling_method=factor(labeling_method, levels=c("metadata", "expression")),
           sex=as.factor(sex)) %>%
    unique() 
  
  ggplot(flow_freq_counts,
         aes(x = labeling_method, 
             stratum = sex, alluvium = row_id,
             y = Freq,
             fill = sex, label = sex)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3) +
    xlab("Label source")+ylab("Number of samples")+
    theme_bw() + theme( panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) + 
    theme(legend.position = "none") 
  
}
alluvPlotSample(human2 %>% select(-pred) %>% rename(gsm=acc))
alluvPlotSample(mouse2 %>% select(-pred) %>% rename(gsm=acc))
alluvPlotSample(rat2 %>% select(-pred) %>% rename(gsm=acc))
## // TODO fix coloring

# study level

alluvPlotStudySumm <- function(ds){
  ggplot(ds,
       aes(x = labeling_method, 
           stratum = sex, 
           alluvium = study,
           y = freq,
           fill = sex, label = sex)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  xlab("Label source")+ylab("Number of studies")+
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 

}
alluvPlotStudySumm(h_summ)
alluvPlotStudySumm(m_summ)
alluvPlotStudySumm(r_summ)

# ------ summarize this all a little better ----- #
h_counts <- h_summ %>% 
  group_by(labeling_method, sex) %>% 
  count() %>%
  ungroup() %>% 
  mutate(organism="human") 
m_counts <- m_summ %>% group_by(labeling_method, sex) %>% count() %>% ungroup() %>% mutate(organism="mouse")
r_counts <- r_summ %>% group_by(labeling_method, sex) %>% count() %>% ungroup() %>% mutate(organism="rat")
all_counts <- do.call(rbind, list(h_counts, m_counts, r_counts))

ggplot(all_counts, aes(x=organism, y=n))+geom_bar(aes(fill=sex), stat="identity")+facet_wrap(.~labeling_method)

all_frac <- do.call(rbind, list(h_counts %>% mutate(n=n/sum(h_counts$n)), 
                                m_counts %>% mutate(n=n/sum(m_counts$n)), 
                                r_counts %>% mutate(n=n/sum(r_counts$n)))) %>%
  mutate(labeling_method=ifelse(as.character(labeling_method)=="exprsex", "expression", "metadata")) %>%
  rename(study_sex=sex)%>%
  mutate(labeling_method=factor(labeling_method, c("metadata", "expression")))

require('viridis')
ggplot(all_frac, aes(x=organism, y=n))+geom_bar(aes(fill=study_sex), stat="identity")+
  facet_wrap(.~labeling_method)+
  ylab("fraction of studies") +
  theme_bw() + theme( panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
  scale_color_viridis(discrete=TRUE, option="D")+
  scale_fill_viridis(discrete=TRUE)
  
ggsave("figures/study_frac_organism.png", dpi="print")
