require('limma')

# TODO: remove rows with zero/low counts
genes <- pos$gene_name
pos_sm <- pos[,2:200]  
neg_sm <- neg[,2:200] 

#dge <- DGEList(counts=counts)
#keep <- filterByExpr(dge, design)
#dge2 <- dge[keep,,]
eset <- cbind(pos_sm,neg_sm)
num_zeros <- apply(eset, 1, function(x) sum(x==0)) # remove > 70% zeros
keep_genes <- (num_zeros/ncol(eset) <= 0.7)
eset2 <- eset[keep_genes,]


design <-as.matrix(c (rep(1, ncol(pos_sm)), rep(0, ncol(neg_sm))))



#dge <- DGEList(counts=counts)
#dge <- calcNormFactors(dge)

rownames(eset2) <- genes[keep_genes]
v <- voom(eset2, design, plot=TRUE) #, normalize="quantile"

fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef=ncol(design))