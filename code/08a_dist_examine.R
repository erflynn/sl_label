require('data.table')
require('tidyverse')
require('ggdendro')
require('gridExtra')

res <- fread("data/cell_dist_mat.txt", data.table=FALSE)
cell_lab <- read_csv("data/cell_line_no_trt.csv")
cell_df <- read_csv("data/00_db_data/cell_syn_df.csv")
# put labels on things
colnames(res) <- cell_lab$acc
rownames(res) <- cell_lab$acc
rownames(cell_lab) <- cell_lab$acc

# pick the shortest description of the cell line that is not numeric
cell_lab2 <- cell_lab %>% left_join(cell_df %>% filter(is.na(as.numeric(cl))) %>%
                                                         group_by(accession) %>% 
                                      top_n(-1, wt=numchar)  %>% select(accession, cl))

# do some hierarchical clustering
# may want to try this melt_dist(res)
samp_id <- sample(1:ncol(res), 1000) %>% sort()
res2 <- res[samp_id,samp_id]
rownames(res2) <- cell_lab$acc[samp_id]
dist_res_sm <- as.dist(res2, upper=TRUE)
cell_lab_sm <- cell_lab2[samp_id,]
res3 <- matrix(0, 1000, 1000)
res3[upper.tri(res3)] <- res2[upper.tri(res2)]
res3[lower.tri(res3)]  <- t(res2)[lower.tri(res2)]
rownames(res3) <- cell_lab$acc[samp_id]
colnames(res3) <- cell_lab$acc[samp_id]

#dist_res <- as.dist(res)
hclust_avg <- hclust(as.dist(res3), method="average")
plot(hclust_avg, labels=FALSE)
ordered_lab <- data.frame("sample_id"=hclust_avg$labels[hclust_avg$order])
ordered_lab2 <- ordered_lab %>% left_join(cell_lab2 %>% select(acc, cl) %>% unique(), by=c("sample_id"="acc"))
ordered_lab3 <- ordered_lab2 %>% mutate(sample_id=factor(sample_id, 
                                                         levels=hclust_avg$labels[hclust_avg$order]), 
                                        cl=as.factor(cl))%>%rename(cell_line=cl)

# sex labels??
human_meta <- read_csv("data/human_metadata2.csv")
head(human_meta)
ordered_lab4 <- ordered_lab3 %>% left_join(human_meta %>% select(acc, sex_lab), by=c("sample_id"="acc"))


require('dendextend')
dend <- hclust_avg %>% as.dendrogram()
col_vec <- c("red", "green", "pink", "blue", "orange", 
             "brown", "yellow", "dark green", "gray", "purple")
cols <- col_vec[ordered_lab3$cell_line]
cols2 <- c("red", "blue")[factor(ordered_lab4$sex_lab)]
plot(dend, leaflab="none")
colored_bars(cbind(cols,cols2), dend, y_shift=-1, rowLabels="", sort_by_labels_order = FALSE)
legend("topright", legend=levels(ordered_lab3$cell_line), col=col_vec, fill=col_vec)
# average distance by cell line

#plot(hclust_avg, labels=FALSE)

# ----- this didnt work ---- #
p1 <- ggdendrogram(hclust_avg, rotate=FALSE)
p2 <- ggplot(ordered_lab3, aes(sample_id, y=1, fill=factor(cell_line)))+geom_tile()+
  scale_y_continuous(expand=c(0,0))+
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")
library(gridExtra)

gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)  

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5))
