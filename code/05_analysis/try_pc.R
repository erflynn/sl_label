# try the PCX and the PCY version
# /scratch/users/erflynn/refine_bio/data/05_train_df/human_rnaseq_sex_neg_expr.RData 
# /scratch/users/erflynn/refine_bio/data/05_train_df/human_rnaseq_sex_neg_expr.RData 

load("data/07_model_dat/human_rnaseq_sex_train_mat.RData")
xy_genes <- read_csv("data/xy_genes_rnaseq.csv")
dim(X_train)
head(xy_genes)
head(X_train[,1:5])
list_genes <- colnames(X_train)
xchr_g <- xy_genes %>% filter(chromosome_name=="X" & transcript %in% list_genes) %>% pull(transcript)
ychr_g <- xy_genes %>% filter(chromosome_name=="Y" & transcript %in% list_genes) %>% pull(transcript)

my_samples <- sample(1:nrow(X_train), 500)


plot(pc_y$x[,"PC1"], pc_x$x[,"PC1"], col=c("red", "blue")[Y_train[my_samples]])


boxplot(t(X_train[my_samples,]))
plot(density(X_train[1,]))

# .... it doesn't help -- is that slightly reassuring? IDK #


X_train1 <- X_train[my_samples,]
Y_train2 <- Y_train[my_samples]
X_valid1 <- X_train[-my_samples,]
Y_valid2 <- Y_train[-my_samples]

require('limma')  
design <- as.matrix(as.numeric(as.character(Y_train2)))
fit <- lmFit(t(X_train1), design)
fit <- eBayes(fit)
tt <-topTable(fit, number=ncol(X_train1))
tt$transcript <- rownames(tt)
x_tt <- tt %>% filter(transcript %in% xchr_g) 
y_tt <- tt %>% filter(transcript %in% ychr_g) 
x_chr_g2 <- head(x_tt$transcript, 10)
y_chr_g2 <- head(y_tt$transcript, 10)

pc_y <- prcomp(X_train[my_samples,y_chr_g2])
pc_x <- prcomp(X_train[my_samples,x_chr_g2])
plot(pc_y$x[,"PC1"], pc_x$x[,"PC1"], 
     col=c("red", "blue")[Y_train[my_samples]])


valid_y <- predict(pc_y, X_valid1[,y_chr_g2])
valid_x <- predict(pc_x, X_valid1[,x_chr_g2])
plot(valid_y[,"PC1"], valid_x[,"PC1"], 
     col=c("red", "blue")[Y_valid2])

table(Y_valid2[valid_y[,"PC1"] < -1]) # 
table(Y_valid2[valid_y[,"PC1"] > -1]) # 
# 86.7% accuracy

# 0. remove replicates of studies

# now do DE analysis

# I'm not sure the normalization is correct!

# 1. start with TMM
#  limma-voom
#  edgeR



# 2. start with counts
# --> try full pipeline
# --> try deseq


# 3. try CNV-kit

# 4. feature selection?

# 5. SVM, RF, LR, PCR again




